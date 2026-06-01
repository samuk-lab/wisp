from __future__ import annotations

import argparse
import gzip
import json
import re
import shutil
import subprocess
from pathlib import Path

import pytest

from sprite_mask import __version__
from sprite_mask.cli import HELP_BANNER, build_parser, main
from sprite_mask.models import WorkflowOutputs

SPRITE_PROGRESS_RE = re.compile(r"^\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2} \[sprite\] Analysis ")


def test_root_help_starts_with_banner_and_version() -> None:
    help_text = build_parser().format_help()

    assert help_text.startswith(f"{HELP_BANNER}\nsprite {__version__}\n\nusage:")


def test_from_alignments_help_describes_fast_mode_replacement() -> None:
    parser = build_parser()
    subparsers = next(
        action for action in parser._actions if isinstance(action, argparse._SubParsersAction)
    )

    help_text = subparsers.choices["from-alignments"].format_help()
    normalized_help = " ".join(help_text.split())

    assert "--fast-mode" in help_text
    assert "strict per-base depth counting is the default" in normalized_help
    assert "--strict-depth" not in help_text


def assert_sprite_progress(log_output: str, message: str) -> None:
    matching_lines = [line for line in log_output.splitlines() if message in line]
    assert matching_lines
    assert SPRITE_PROGRESS_RE.match(matching_lines[0])


def test_main_all_sites_vcf_writes_indexed_population_bed(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    missing = [tool for tool in ("bgzip", "tabix") if shutil.which(tool) is None]
    if missing:
        pytest.skip(f"external VCF workflow tool(s) unavailable: {', '.join(missing)}")

    vcf = tmp_path / "all_sites.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n"
        "chr1\t1\t.\tA\t.\t.\t.\t.\tGT:DP\t0/0:5\t0/0:2\n"
        "chr1\t2\t.\tC\t.\t.\t.\t.\tGT:DP\t0/0:6\t0/0:7\n"
        "chr1\t3\t.\tG\t.\t.\t.\t.\tGT:DP\t0/0:8\t0/0:9\n"
        "chr1\t4\t.\tT\t.\t.\t.\t.\tGT:DP\t0/0:0\t0/0:0\n"
        "chr1\t5\t.\tA\t.\t.\t.\t.\tGT:DP\t0/0:5\t0/0:5\n"
    )
    popfile = tmp_path / "popfile.tsv"
    popfile.write_text("sample_id\tpopulation\ns1\tpopA\ns2\tpopB\n")
    out_dir = tmp_path / "out"
    work_dir = tmp_path / "work"

    status = main(
        [
            "from-vcf",
            "--all-sites-vcf",
            str(vcf),
            "--popfile",
            str(popfile),
            "--min-dp",
            "5",
            "--out",
            str(out_dir),
            "--work",
            str(work_dir),
        ]
    )

    assert status == 0
    captured = capsys.readouterr()
    assert captured.out == ""
    assert_sprite_progress(captured.err, "Analysis start: validating VCF workflow inputs")
    assert_sprite_progress(captured.err, "Analysis VCF: building population counts")
    assert_sprite_progress(captured.err, "Analysis complete: wrote")
    assert "sprite.bed.gz" in captured.err

    population_bed = out_dir / "sprite.bed.gz"
    population_index = Path(f"{population_bed}.tbi")
    assert population_bed.exists()
    assert population_index.exists()
    assert not work_dir.exists()

    with gzip.open(population_bed, "rt") as handle:
        lines = handle.read().splitlines()
    metadata = json.loads(lines[0].split("\t", maxsplit=1)[1])
    assert metadata["source_mode"] == "all_sites_vcf"
    assert metadata["columns"] == ["chrom", "start", "end", "popA", "popB"]
    assert metadata["population_sample_counts"] == {"popA": 1, "popB": 1}
    assert lines[1:] == [
        "#chrom\tstart\tend\tpopA\tpopB",
        "chr1\t0\t1\t1\t0",
        "chr1\t1\t3\t1\t1",
        "chr1\t4\t5\t1\t1",
    ]

    completed = subprocess.run(
        ["tabix", str(population_bed), "chr1:2-3"],
        check=True,
        capture_output=True,
        text=True,
    )
    assert completed.stdout.strip() == "chr1\t1\t3\t1\t1"


def test_main_builds_alignment_run_config(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    seen_config = None

    def fake_run_workflow(config: object) -> WorkflowOutputs:
        nonlocal seen_config
        seen_config = config
        return WorkflowOutputs(
            population_count_bed_gz=tmp_path / "out" / "sprite.bed.gz",
            population_count_bed_index=tmp_path / "out" / "sprite.bed.gz.tbi",
        )

    monkeypatch.setattr("sprite_mask.cli.run_workflow", fake_run_workflow)

    status = main(
        [
            "from-alignments",
            "--samples",
            str(tmp_path / "samples.tsv"),
            "--min-dp",
            "30",
            "--out",
            str(tmp_path / "out"),
            "--work",
            str(tmp_path / "work"),
            "--threads",
            "4",
            "--jobs",
            "2",
            "--mask",
            str(tmp_path / "targets.bed"),
            "--variants-vcf",
            str(tmp_path / "variants.vcf.gz"),
            "--min-mapq",
            "20",
            "--exclude-flag",
            "1796",
            "--reference",
            str(tmp_path / "ref.fa"),
            "--fast-mode",
            "--keep-work",
            "--force",
        ]
    )

    assert status == 0
    assert seen_config is not None
    assert seen_config.samples_path == tmp_path / "samples.tsv"
    assert seen_config.min_dp == 30
    assert seen_config.output_prefix == "sprite"
    assert seen_config.threads == 4
    assert seen_config.jobs == 2
    assert seen_config.mask_bed == tmp_path / "targets.bed"
    assert seen_config.variants_vcf == tmp_path / "variants.vcf.gz"
    assert seen_config.min_mapq == 20
    assert seen_config.exclude_flag == 1796
    assert seen_config.reference == tmp_path / "ref.fa"
    assert seen_config.fast_mode is True
    assert seen_config.keep_work is True
    assert seen_config.force is True

    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == ""


def test_main_builds_alignment_run_config_with_variants_vcf_default_min_dp(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    seen_config = None

    def fake_run_workflow(config: object) -> WorkflowOutputs:
        nonlocal seen_config
        seen_config = config
        return WorkflowOutputs(
            population_count_bed_gz=tmp_path / "out" / "sprite.bed.gz",
            population_count_bed_index=tmp_path / "out" / "sprite.bed.gz.tbi",
        )

    monkeypatch.setattr("sprite_mask.cli.run_workflow", fake_run_workflow)

    status = main(
        [
            "from-alignments",
            "--samples",
            str(tmp_path / "samples.tsv"),
            "--variants-vcf",
            str(tmp_path / "variants.vcf.gz"),
            "--out",
            str(tmp_path / "out"),
        ]
    )

    assert status == 0
    assert seen_config is not None
    assert seen_config.min_dp is None
    assert seen_config.variants_vcf == tmp_path / "variants.vcf.gz"


def test_main_accepts_output_prefix(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    seen_config = None

    def fake_run_workflow(config: object) -> WorkflowOutputs:
        nonlocal seen_config
        seen_config = config
        return WorkflowOutputs(
            population_count_bed_gz=tmp_path / "out" / "custom.bed.gz",
            population_count_bed_index=tmp_path / "out" / "custom.bed.gz.tbi",
        )

    monkeypatch.setattr("sprite_mask.cli.run_workflow", fake_run_workflow)

    status = main(
        [
            "from-vcf",
            "--all-sites-vcf",
            str(tmp_path / "all_sites.vcf.gz"),
            "--popfile",
            str(tmp_path / "popfile.tsv"),
            "--min-dp",
            "30",
            "--out",
            str(tmp_path / "out"),
            "--output-prefix",
            "custom",
        ]
    )

    assert status == 0
    assert seen_config is not None
    assert seen_config.output_prefix == "custom"


def test_main_from_vcf_accepts_snps_only(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    seen_config = None

    def fake_run_workflow(config: object) -> WorkflowOutputs:
        nonlocal seen_config
        seen_config = config
        return WorkflowOutputs(
            population_count_bed_gz=tmp_path / "out" / "sprite.bed.gz",
            population_count_bed_index=tmp_path / "out" / "sprite.bed.gz.tbi",
        )

    monkeypatch.setattr("sprite_mask.cli.run_workflow", fake_run_workflow)

    status = main(
        [
            "from-vcf",
            "--all-sites-vcf",
            str(tmp_path / "all_sites.vcf.gz"),
            "--popfile",
            str(tmp_path / "popfile.tsv"),
            "--min-dp",
            "30",
            "--out",
            str(tmp_path / "out"),
            "--snps-only",
        ]
    )

    assert status == 0
    assert seen_config is not None
    assert seen_config.snps_only is True


def test_main_reports_subprocess_errors(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def fake_run_workflow(_config: object) -> WorkflowOutputs:
        raise subprocess.CalledProcessError(9, ["tool", "arg"], stderr="bad things\n")

    monkeypatch.setattr("sprite_mask.cli.run_workflow", fake_run_workflow)

    status = main(
        [
            "from-alignments",
            "--samples",
            str(tmp_path / "samples.tsv"),
            "--min-dp",
            "30",
            "--out",
            str(tmp_path / "out"),
        ]
    )

    assert status == 1
    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == (
        "sprite: command failed with exit code 9\n"
        "tool arg\n"
        "bad things\n"
    )


def test_main_reports_regular_exceptions(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def fake_run_workflow(_config: object) -> WorkflowOutputs:
        raise ValueError("bad config")

    monkeypatch.setattr("sprite_mask.cli.run_workflow", fake_run_workflow)

    status = main(
        [
            "from-alignments",
            "--samples",
            str(tmp_path / "samples.tsv"),
            "--min-dp",
            "30",
            "--out",
            str(tmp_path / "out"),
        ]
    )

    assert status == 1
    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == "sprite: error: bad config\n"


def test_main_no_subcommand_returns_error(capsys: pytest.CaptureFixture[str]) -> None:
    status = main([])

    assert status == 1
    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err.startswith(f"{HELP_BANNER}\nsprite {__version__}\n\nusage:")


def test_main_dry_run_skips_execution_and_returns_zero(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls: list[object] = []

    def fake_run_workflow(config: object) -> WorkflowOutputs:
        calls.append(config)
        return WorkflowOutputs(
            population_count_bed_gz=tmp_path / "out" / "sprite.bed.gz",
            population_count_bed_index=tmp_path / "out" / "sprite.bed.gz.tbi",
        )

    monkeypatch.setattr("sprite_mask.cli.run_workflow", fake_run_workflow)

    status = main(
        [
            "from-alignments",
            "--samples",
            str(tmp_path / "samples.tsv"),
            "--min-dp",
            "30",
            "--out",
            str(tmp_path / "out"),
            "--dry-run",
        ]
    )

    assert status == 0
    assert len(calls) == 1
    assert calls[0].dry_run is True  # type: ignore[union-attr]


def test_main_debug_reraises_exception(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def fake_run_workflow(_config: object) -> WorkflowOutputs:
        raise ValueError("internal error")

    monkeypatch.setattr("sprite_mask.cli.run_workflow", fake_run_workflow)

    with pytest.raises(ValueError, match="internal error"):
        main(
            [
                "from-alignments",
                "--samples",
                str(tmp_path / "samples.tsv"),
                "--min-dp",
                "30",
                "--out",
                str(tmp_path / "out"),
                "--debug",
            ]
        )
