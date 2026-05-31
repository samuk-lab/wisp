from __future__ import annotations

import gzip
import subprocess
from collections.abc import Sequence
from pathlib import Path
from typing import TextIO

import pytest

from sprite_mask.config import AlignmentRunConfig, VcfRunConfig
from sprite_mask.models import MosdepthOutputs, Sample
from sprite_mask.validation import validate_vcf_inputs
from sprite_mask.workflow import (
    _build_from_alignments,
    _build_from_all_sites_vcf,
    _cleanup_work_files,
    _make_sample_pass_bed,
    _make_sample_pass_beds,
    _prepare_targets,
    _required_tools,
    _sort_bgzip_tabix_bed,
    run_workflow,
)


def test_validate_vcf_inputs_rejects_missing_vcf(tmp_path: Path) -> None:
    popfile = tmp_path / "popfile.tsv"
    popfile.write_text("sample_id\tpopulation\n")

    with pytest.raises(ValueError, match="does not exist"):
        validate_vcf_inputs(tmp_path / "missing.vcf", popfile)


def test_validate_vcf_inputs_rejects_vcf_directory(tmp_path: Path) -> None:
    vcf_dir = tmp_path / "vcf_dir"
    vcf_dir.mkdir()
    popfile = tmp_path / "popfile.tsv"
    popfile.write_text("sample_id\tpopulation\n")

    with pytest.raises(ValueError, match="is a directory"):
        validate_vcf_inputs(vcf_dir, popfile)


def test_validate_vcf_inputs_rejects_missing_popfile(tmp_path: Path) -> None:
    vcf = tmp_path / "all_sites.vcf"
    vcf.write_text("vcf")

    with pytest.raises(ValueError, match="--popfile does not exist"):
        validate_vcf_inputs(vcf, tmp_path / "missing.tsv")


def test_validate_vcf_inputs_rejects_popfile_directory(tmp_path: Path) -> None:
    vcf = tmp_path / "all_sites.vcf"
    vcf.write_text("vcf")
    popfile_dir = tmp_path / "popfile_dir"
    popfile_dir.mkdir()

    with pytest.raises(ValueError, match="--popfile is a directory"):
        validate_vcf_inputs(vcf, popfile_dir)


def test_required_tools_for_vcf_mode_omits_mosdepth_and_bedtools(tmp_path: Path) -> None:
    vcf = tmp_path / "all_sites.vcf"
    vcf.write_text("vcf")
    config = VcfRunConfig(
        all_sites_vcf=vcf,
        popfile_path=tmp_path / "popfile.tsv",
        min_dp=5,
        out_dir=tmp_path / "out",
    )
    assert _required_tools(config) == ("bgzip", "tabix")


def test_required_tools_for_threshold_zero_omits_mosdepth(tmp_path: Path) -> None:
    config = AlignmentRunConfig(
        samples_path=tmp_path / "samples.tsv",
        min_dp=0,
        out_dir=tmp_path / "out",
        mask_bed=tmp_path / "targets.bed",
    )

    assert _required_tools(config) == ("samtools", "bedtools", "bgzip", "tabix")


def test_required_tools_for_alignment_mode_includes_mosdepth(tmp_path: Path) -> None:
    config = AlignmentRunConfig(
        samples_path=tmp_path / "samples.tsv",
        min_dp=10,
        out_dir=tmp_path / "out",
    )

    assert _required_tools(config) == ("samtools", "mosdepth", "bedtools", "bgzip", "tabix")


def test_run_workflow_alignment_mode_dispatches_and_cleans_work_files(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    sample = Sample("s1", "popA", tmp_path / "s1.bam")
    config = AlignmentRunConfig(
        samples_path=tmp_path / "samples.tsv",
        min_dp=10,
        out_dir=tmp_path / "out",
        work_dir=tmp_path / "work",
    )
    calls: list[tuple[list[Sample], AlignmentRunConfig]] = []

    monkeypatch.setattr("sprite_mask.workflow.read_samples", lambda _path: [sample])
    monkeypatch.setattr("sprite_mask.workflow.require_executables", lambda _names: None)
    monkeypatch.setattr(
        "sprite_mask.workflow.validate_alignment_sample_headers",
        lambda _samples: None,
    )
    monkeypatch.setattr(
        "sprite_mask.workflow._build_from_all_sites_vcf",
        lambda *_args: pytest.fail("VCF workflow should not be called"),
    )

    def fake_build_from_alignments(
        samples: list[Sample],
        config_arg: AlignmentRunConfig,
        generated_work_files: list[Path],
    ) -> None:
        calls.append((samples, config_arg))
        generated = config_arg.resolved_work_dir / "generated.tmp"
        generated.write_text("generated")
        generated_work_files.append(generated)

    monkeypatch.setattr(
        "sprite_mask.workflow._build_from_alignments",
        fake_build_from_alignments,
    )

    outputs = run_workflow(config)

    assert outputs.population_count_bed_gz == tmp_path / "out" / "sprite.bed.gz"
    assert calls == [([sample], config)]
    assert not config.resolved_work_dir.exists()


def test_prepare_targets_normalizes_sorts_and_tracks_outputs(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls: list[tuple[str, Path, Path]] = []

    def fake_normalize(targets_bed: Path, out_bed: Path) -> Path:
        calls.append(("normalize", targets_bed, out_bed))
        out_bed.write_text("chr1\t0\t10\n")
        return out_bed

    def fake_sort(in_bed: Path, out_bed: Path) -> Path:
        calls.append(("sort", in_bed, out_bed))
        out_bed.write_text(in_bed.read_text())
        return out_bed

    monkeypatch.setattr("sprite_mask.workflow.normalize_targets_bed", fake_normalize)
    monkeypatch.setattr("sprite_mask.workflow.sort_and_merge_bed", fake_sort)

    config = AlignmentRunConfig(
        samples_path=tmp_path / "samples.tsv",
        min_dp=10,
        out_dir=tmp_path / "out",
        work_dir=tmp_path / "work",
        mask_bed=tmp_path / "targets.bed",
    )
    config.resolved_work_dir.mkdir(parents=True)
    generated: list[Path] = []

    target_bed = _prepare_targets(config, generated)

    normalized = tmp_path / "work" / "targets.3col.bed"
    sorted_merged = tmp_path / "work" / "targets.3col.sorted.merged.bed"
    assert target_bed == sorted_merged
    assert generated == [normalized, sorted_merged]
    assert calls == [
        ("normalize", tmp_path / "targets.bed", normalized),
        ("sort", normalized, sorted_merged),
    ]


def test_prepare_targets_without_targets_returns_none(tmp_path: Path) -> None:
    generated: list[Path] = []

    assert (
        _prepare_targets(
            AlignmentRunConfig(
                samples_path=tmp_path / "samples.tsv",
                min_dp=10,
                out_dir=tmp_path / "out",
            ),
            generated,
        )
        is None
    )
    assert generated == []


def test_make_sample_pass_bed_threshold_zero_copies_target_bed(tmp_path: Path) -> None:
    target_bed = tmp_path / "targets.bed"
    target_bed.write_text("chr1\t0\t10\n")
    config = AlignmentRunConfig(
        samples_path=tmp_path / "samples.tsv",
        min_dp=0,
        out_dir=tmp_path / "out",
        work_dir=tmp_path / "work",
        mask_bed=target_bed,
    )
    config.resolved_work_dir.mkdir(parents=True)

    pass_bed, mosdepth_outputs, generated = _make_sample_pass_bed(
        Sample("s1", "popA", tmp_path / "s1.bam"),
        config,
        target_bed,
    )

    assert mosdepth_outputs is None
    assert pass_bed == tmp_path / "work" / "s1.d0.pass.targets.bed"
    assert pass_bed.read_text() == "chr1\t0\t10\n"
    assert generated == [pass_bed]


def test_make_sample_pass_bed_threshold_zero_requires_targets(tmp_path: Path) -> None:
    config = AlignmentRunConfig(
        samples_path=tmp_path / "samples.tsv",
        min_dp=0,
        out_dir=tmp_path / "out",
        work_dir=tmp_path / "work",
    )

    with pytest.raises(ValueError, match="requires --mask"):
        _make_sample_pass_bed(Sample("s1", "popA", tmp_path / "s1.bam"), config, None)


def test_make_sample_pass_bed_without_targets_uses_merged_pass_bed(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    config = AlignmentRunConfig(
        samples_path=tmp_path / "samples.tsv",
        min_dp=10,
        out_dir=tmp_path / "out",
        work_dir=tmp_path / "work",
    )
    config.resolved_work_dir.mkdir(parents=True)
    outputs = _mosdepth_outputs(tmp_path / "work" / "s1.d10")

    def fake_run_mosdepth(_sample: Sample, _config: AlignmentRunConfig) -> MosdepthOutputs:
        return outputs

    def fake_extract(quantized_bed_gz: Path, out_bed: Path) -> Path:
        assert quantized_bed_gz == outputs.quantized_bed_gz
        out_bed.write_text("chr1\t0\t10\n")
        return out_bed

    monkeypatch.setattr("sprite_mask.workflow.run_mosdepth", fake_run_mosdepth)
    monkeypatch.setattr("sprite_mask.workflow.extract_merged_pass_intervals", fake_extract)

    pass_bed, mosdepth_outputs, generated = _make_sample_pass_bed(
        Sample("s1", "popA", tmp_path / "s1.bam"),
        config,
        None,
    )

    assert pass_bed == tmp_path / "work" / "s1.d10.pass.bed"
    assert pass_bed.read_text() == "chr1\t0\t10\n"
    assert mosdepth_outputs == outputs
    assert generated == [
        outputs.quantized_bed_gz,
        outputs.quantized_bed_index,
        outputs.summary,
        outputs.global_dist,
        outputs.stderr_log,
        pass_bed,
    ]


def test_make_sample_pass_bed_with_targets_clips_merged_pass_bed(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    config = AlignmentRunConfig(
        samples_path=tmp_path / "samples.tsv",
        min_dp=10,
        out_dir=tmp_path / "out",
        work_dir=tmp_path / "work",
    )
    config.resolved_work_dir.mkdir(parents=True)
    target_bed = tmp_path / "targets.bed"
    outputs = _mosdepth_outputs(tmp_path / "work" / "s1.d10")
    calls: list[tuple[Path, Path, Path]] = []

    monkeypatch.setattr(
        "sprite_mask.workflow.run_mosdepth",
        lambda _sample, _config: outputs,
    )
    monkeypatch.setattr(
        "sprite_mask.workflow.extract_merged_pass_intervals",
        lambda _quantized, out_bed: out_bed.write_text("chr1\t0\t20\n") or out_bed,
    )

    def fake_intersect(merged_pass_bed: Path, targets_bed: Path, out_bed: Path) -> Path:
        calls.append((merged_pass_bed, targets_bed, out_bed))
        out_bed.write_text("chr1\t5\t10\n")
        return out_bed

    monkeypatch.setattr("sprite_mask.workflow.intersect_sort_merge", fake_intersect)

    pass_bed, returned_outputs, generated = _make_sample_pass_bed(
        Sample("s1", "popA", tmp_path / "s1.bam"),
        config,
        target_bed,
    )

    merged_pass_bed = tmp_path / "work" / "s1.d10.pass.bed"
    clipped_pass_bed = tmp_path / "work" / "s1.d10.pass.targets.bed"
    assert pass_bed == clipped_pass_bed
    assert returned_outputs == outputs
    assert pass_bed.read_text() == "chr1\t5\t10\n"
    assert calls == [(merged_pass_bed, target_bed, clipped_pass_bed)]
    assert generated[-2:] == [merged_pass_bed, clipped_pass_bed]


def test_make_sample_pass_beds_uses_parallel_jobs_and_tracks_work_files(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    samples = [
        Sample("s1", "popA", tmp_path / "s1.bam"),
        Sample("s2", "popB", tmp_path / "s2.bam"),
    ]
    config = AlignmentRunConfig(
        samples_path=tmp_path / "samples.tsv",
        min_dp=10,
        out_dir=tmp_path / "out",
        work_dir=tmp_path / "work",
        jobs=4,
    )

    def fake_make_sample_pass_bed(
        sample: Sample,
        _config: AlignmentRunConfig,
        _target_bed: Path | None,
    ) -> tuple[Path, None, list[Path]]:
        pass_bed = tmp_path / f"{sample.sample_id}.pass.bed"
        sample_log = tmp_path / f"{sample.sample_id}.log"
        return pass_bed, None, [sample_log]

    monkeypatch.setattr("sprite_mask.workflow._make_sample_pass_bed", fake_make_sample_pass_bed)
    generated: list[Path] = []

    pass_beds = _make_sample_pass_beds(samples, config, None, generated)

    assert pass_beds == [tmp_path / "s1.pass.bed", tmp_path / "s2.pass.bed"]
    assert generated == [tmp_path / "s1.log", tmp_path / "s2.log"]


def test_make_sample_pass_beds_uses_sequential_path_for_one_job(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    samples = [
        Sample("s1", "popA", tmp_path / "s1.bam"),
        Sample("s2", "popB", tmp_path / "s2.bam"),
    ]
    config = AlignmentRunConfig(
        samples_path=tmp_path / "samples.tsv",
        min_dp=10,
        out_dir=tmp_path / "out",
        work_dir=tmp_path / "work",
        jobs=1,
    )
    visited: list[str] = []

    def fake_make_sample_pass_bed(
        sample: Sample,
        _config: AlignmentRunConfig,
        _target_bed: Path | None,
    ) -> tuple[Path, None, list[Path]]:
        visited.append(sample.sample_id)
        return tmp_path / f"{sample.sample_id}.pass.bed", None, []

    monkeypatch.setattr("sprite_mask.workflow._make_sample_pass_bed", fake_make_sample_pass_bed)

    pass_beds = _make_sample_pass_beds(samples, config, None, [])

    assert visited == ["s1", "s2"]
    assert pass_beds == [tmp_path / "s1.pass.bed", tmp_path / "s2.pass.bed"]


def test_build_from_alignments_uses_single_input_multiinter_for_one_sample(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    config = AlignmentRunConfig(
        samples_path=tmp_path / "samples.tsv",
        min_dp=10,
        out_dir=tmp_path / "out",
        work_dir=tmp_path / "work",
    )
    config.resolved_work_dir.mkdir(parents=True)
    sample = Sample("s1", "popA", tmp_path / "s1.bam")
    pass_bed = tmp_path / "s1.pass.bed"
    calls: list[str] = []

    monkeypatch.setattr(
        "sprite_mask.workflow._prepare_targets",
        lambda _config, _generated: None,
    )
    monkeypatch.setattr(
        "sprite_mask.workflow._make_sample_pass_beds",
        lambda _samples, _config, _target_bed, _generated: [pass_bed],
    )

    def fake_single(pass_bed_arg: Path, name: str, out_tsv: Path) -> Path:
        calls.append("single")
        assert pass_bed_arg == pass_bed
        assert name == "s1"
        out_tsv.write_text("chrom\tstart\tend\tnum\tlist\ts1\n")
        return out_tsv

    def fake_collapse(
        samples: list[Sample],
        multiinter_tsv: Path,
        output_bed: Path,
        *,
        metadata: dict[str, object],
    ) -> Path:
        calls.append("collapse")
        assert samples == [sample]
        assert multiinter_tsv.name == "cohort.d10.multiinter.tsv"
        assert metadata["samples_path"] == str(config.samples_path)
        output_bed.write_text("#chrom\tstart\tend\tpopA\n")
        return output_bed

    def fake_sort(in_bed: Path, out_bed_gz: Path) -> None:
        calls.append("sort")
        assert in_bed.name == "cohort.d10.population_count_quantized.bed"
        assert out_bed_gz == tmp_path / "out" / "sprite.bed.gz"

    monkeypatch.setattr("sprite_mask.workflow.write_single_input_multiinter", fake_single)
    monkeypatch.setattr(
        "sprite_mask.workflow.run_multiinter",
        lambda *_args: pytest.fail("run_multiinter should not be called"),
    )
    monkeypatch.setattr("sprite_mask.workflow.collapse_population_counts", fake_collapse)
    monkeypatch.setattr("sprite_mask.workflow._sort_bgzip_tabix_bed", fake_sort)
    generated: list[Path] = []

    _build_from_alignments([sample], config, generated)

    assert calls == ["single", "collapse", "sort"]
    assert generated == [
        tmp_path / "work" / "cohort.d10.multiinter.tsv",
        tmp_path / "work" / "cohort.d10.population_count_quantized.bed",
    ]


def test_build_from_alignments_uses_bedtools_multiinter_for_multiple_samples(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    config = AlignmentRunConfig(
        samples_path=tmp_path / "samples.tsv",
        min_dp=10,
        out_dir=tmp_path / "out",
        work_dir=tmp_path / "work",
    )
    config.resolved_work_dir.mkdir(parents=True)
    samples = [
        Sample("s1", "popA", tmp_path / "s1.bam"),
        Sample("s2", "popB", tmp_path / "s2.bam"),
    ]
    pass_beds = [tmp_path / "s1.pass.bed", tmp_path / "s2.pass.bed"]
    calls: list[tuple[list[Path], list[str], Path]] = []

    monkeypatch.setattr(
        "sprite_mask.workflow._prepare_targets",
        lambda _config, _generated: None,
    )
    monkeypatch.setattr(
        "sprite_mask.workflow._make_sample_pass_beds",
        lambda _samples, _config, _target_bed, _generated: pass_beds,
    )

    def fake_multi(pass_beds_arg: Sequence[Path], names: Sequence[str], out_tsv: Path) -> Path:
        calls.append((list(pass_beds_arg), list(names), out_tsv))
        out_tsv.write_text("chrom\tstart\tend\tnum\tlist\ts1\ts2\n")
        return out_tsv

    monkeypatch.setattr("sprite_mask.workflow.run_multiinter", fake_multi)
    monkeypatch.setattr(
        "sprite_mask.workflow.collapse_population_counts",
        lambda _samples, _multiinter, output_bed, *, metadata: output_bed.write_text(
            "#chrom\tstart\tend\tpopA\tpopB\n"
        )
        or output_bed,
    )
    monkeypatch.setattr("sprite_mask.workflow._sort_bgzip_tabix_bed", lambda *_args: None)

    _build_from_alignments(samples, config, [])

    assert calls == [
        (
            pass_beds,
            ["s1", "s2"],
            tmp_path / "work" / "cohort.d10.multiinter.tsv",
        )
    ]


def test_build_from_all_sites_vcf_writes_population_bed_and_metadata(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    vcf = tmp_path / "all_sites.vcf"
    vcf.write_text("vcf")
    popfile = tmp_path / "popfile.tsv"
    popfile.write_text("sample_id\tpopulation\n")
    targets = tmp_path / "targets.bed"
    samples = [Sample("s1", "popA")]
    config = VcfRunConfig(
        all_sites_vcf=vcf,
        popfile_path=popfile,
        min_dp=5,
        out_dir=tmp_path / "out",
        work_dir=tmp_path / "work",
        mask_bed=targets,
    )
    config.resolved_work_dir.mkdir(parents=True)
    calls: dict[str, object] = {}

    def fake_build(
        samples_arg: list[Sample],
        all_sites_vcf: Path,
        output_bed: Path,
        *,
        threshold: int,
        targets_bed: Path | None,
        snps_only: bool,
        metadata: dict[str, object],
    ) -> Path:
        calls["build"] = (
            samples_arg,
            all_sites_vcf,
            output_bed,
            threshold,
            targets_bed,
            snps_only,
            metadata,
        )
        output_bed.write_text("#chrom\tstart\tend\tpopA\n")
        return output_bed

    def fake_sort(in_bed: Path, out_bed_gz: Path) -> None:
        calls["sort"] = (in_bed, out_bed_gz)

    monkeypatch.setattr(
        "sprite_mask.workflow.build_population_counts_from_all_sites_vcf",
        fake_build,
    )
    monkeypatch.setattr("sprite_mask.workflow._sort_bgzip_tabix_bed", fake_sort)
    generated: list[Path] = []

    _build_from_all_sites_vcf(samples, config, generated)

    output_bed = tmp_path / "work" / "cohort.d5.population_count_quantized.bed"
    assert generated == [output_bed]
    build_call = calls["build"]
    assert build_call[:5] == (samples, vcf, output_bed, 5, targets)
    assert build_call[5] is False
    assert build_call[6] == {
        "min_dp": 5,
        "sample_count": 1,
        "popfile": str(popfile),
        "all_sites_vcf": str(vcf),
        "mask_bed": str(targets),
        "snps_only": False,
    }
    assert calls["sort"] == (output_bed, tmp_path / "out" / "sprite.bed.gz")


def test_sort_bgzip_tabix_bed_preserves_headers_sorts_body_and_removes_temps(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    in_bed = tmp_path / "population_count.bed"
    in_bed.write_text(
        "#sprite_mask_metadata\t{}\n"
        "chrom\tstart\tend\tpopA\n"
        "chr2\t5\t6\t1\n"
        "chr1\t2\t3\t1\n"
        "\n"
        "chr1\t0\t1\t1\n"
    )
    out_bed_gz = tmp_path / "out" / "sprite.bed.gz"

    def fake_run(
        command: list[str],
        *,
        check: bool,
        stdout: TextIO | None = None,
        text: bool | None = None,
        env: dict[str, str] | None = None,
    ) -> subprocess.CompletedProcess[str]:
        assert check is True
        if command[0] == "sort":
            assert stdout is not None
            assert env is not None
            assert env.get("LC_ALL") == "C"
            body_path = Path(command[-1])
            sorted_lines = sorted(
                body_path.read_text().splitlines(),
                key=lambda line: (line.split("\t")[0], int(line.split("\t")[1])),
            )
            stdout.write("\n".join(sorted_lines) + "\n")
            return subprocess.CompletedProcess(command, 0)
        if command[0] == "bgzip":
            sorted_bed = Path(command[-1])
            with gzip.open(out_bed_gz, "wt") as handle:
                handle.write(sorted_bed.read_text())
            return subprocess.CompletedProcess(command, 0)
        if command[0] == "tabix":
            Path(f"{command[-1]}.tbi").write_text("index")
            return subprocess.CompletedProcess(command, 0)
        raise AssertionError(f"unexpected command: {command}")

    monkeypatch.setattr("sprite_mask.workflow.subprocess.run", fake_run)

    _sort_bgzip_tabix_bed(in_bed, out_bed_gz)

    with gzip.open(out_bed_gz, "rt") as handle:
        assert handle.read() == (
            "#sprite_mask_metadata\t{}\n"
            "#chrom\tstart\tend\tpopA\n"
            "chr1\t0\t1\t1\n"
            "chr1\t2\t3\t1\n"
            "chr2\t5\t6\t1\n"
        )
    assert Path(f"{out_bed_gz}.tbi").read_text() == "index"
    assert not (tmp_path / "out" / "sprite.bed").exists()
    assert not (tmp_path / "out" / "sprite.bed.body").exists()
    assert not (tmp_path / "out" / "sprite.bed.sorted_body").exists()


def test_cleanup_work_files_removes_known_files_and_empty_work_dir(tmp_path: Path) -> None:
    work_dir = tmp_path / "work"
    work_dir.mkdir()
    first = work_dir / "first.tmp"
    second = work_dir / "second.tmp"
    first.write_text("first")
    second.write_text("second")

    _cleanup_work_files([first, work_dir / "missing.tmp", second], work_dir)

    assert not work_dir.exists()


def test_cleanup_work_files_leaves_non_empty_work_dir(tmp_path: Path) -> None:
    work_dir = tmp_path / "work"
    work_dir.mkdir()
    generated = work_dir / "generated.tmp"
    survivor = work_dir / "survivor.tmp"
    generated.write_text("generated")
    survivor.write_text("survivor")

    _cleanup_work_files([generated], work_dir)

    assert not generated.exists()
    assert survivor.exists()
    assert work_dir.exists()


def test_full_all_sites_run_workflow_rejects_existing_outputs(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    population_bed = tmp_path / "out" / "sprite.bed.gz"
    population_bed.parent.mkdir()
    population_bed.write_text("old")
    vcf = tmp_path / "all_sites.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n"
    )
    popfile = tmp_path / "popfile.tsv"
    popfile.write_text("sample_id\tpopulation\ns1\tpopA\n")

    monkeypatch.setattr("sprite_mask.workflow.require_executables", lambda _names: None)

    with pytest.raises(FileExistsError, match="pass --force"):
        run_workflow(
            VcfRunConfig(
                all_sites_vcf=vcf,
                popfile_path=popfile,
                min_dp=5,
                out_dir=tmp_path / "out",
            )
        )


def _mosdepth_outputs(prefix: Path) -> MosdepthOutputs:
    return MosdepthOutputs(
        prefix=prefix,
        quantized_bed_gz=Path(f"{prefix}.quantized.bed.gz"),
        quantized_bed_index=Path(f"{prefix}.quantized.bed.gz.csi"),
        summary=Path(f"{prefix}.mosdepth.summary.txt"),
        global_dist=Path(f"{prefix}.mosdepth.global.dist.txt"),
        stderr_log=Path(f"{prefix}.mosdepth.stderr.log"),
    )
