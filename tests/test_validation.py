from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from sprite_mask.models import Sample
from sprite_mask.validation import (
    ensure_parent_dirs,
    refuse_existing_outputs,
    require_executables,
    validate_alignment_sample_headers,
    validate_jobs,
    validate_threads,
    validate_threshold,
    validate_variants_vcf_input,
)


def test_validate_threshold_rejects_negative() -> None:
    with pytest.raises(ValueError, match="non-negative"):
        validate_threshold(-1)


def test_validate_threshold_zero_requires_targets_or_vcf(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="requires --mask"):
        validate_threshold(0)

    validate_threshold(0, targets_bed=tmp_path / "targets.bed")
    validate_threshold(0, all_sites_vcf=tmp_path / "all_sites.vcf")


def test_validate_threads_and_jobs_reject_values_below_one() -> None:
    with pytest.raises(ValueError, match="--threads"):
        validate_threads(0)
    with pytest.raises(ValueError, match="--jobs"):
        validate_jobs(0)

    validate_threads(1)
    validate_jobs(1)


def test_validate_variants_vcf_input_rejects_missing_or_directory(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="--variants-vcf does not exist"):
        validate_variants_vcf_input(tmp_path / "missing.vcf")

    vcf_dir = tmp_path / "variants_dir"
    vcf_dir.mkdir()
    with pytest.raises(ValueError, match="--variants-vcf is a directory"):
        validate_variants_vcf_input(vcf_dir)

    variants_vcf = tmp_path / "variants.vcf"
    variants_vcf.write_text("vcf")
    validate_variants_vcf_input(variants_vcf)


def test_validate_alignment_sample_headers_accepts_matching_read_group_samples(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    alignment = tmp_path / "s1.bam"
    alignment.write_bytes(b"bam")

    def fake_run(
        command: list[str],
        *,
        check: bool,
        capture_output: bool,
        text: bool,
    ) -> object:
        assert command == ["samtools", "view", "-H", str(alignment)]
        assert check is False
        assert capture_output is True
        assert text is True
        return subprocess.CompletedProcess(
            command,
            0,
            stdout="@HD\tVN:1.6\n@RG\tID:rg1\tSM:s1\n@RG\tID:rg2\tSM:s1\n",
            stderr="",
        )

    monkeypatch.setattr("sprite_mask.validation.subprocess.run", fake_run)

    validate_alignment_sample_headers([Sample("s1", "popA", alignment)])


def test_validate_alignment_sample_headers_rejects_mismatched_read_group_samples(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    alignment = tmp_path / "s1.bam"
    alignment.write_bytes(b"bam")

    def fake_run(
        command: list[str],
        *,
        check: bool,
        capture_output: bool,
        text: bool,
    ) -> object:
        return subprocess.CompletedProcess(
            command,
            0,
            stdout="@RG\tID:rg1\tSM:s1\n@RG\tID:rg2\tSM:s2\n",
            stderr="",
        )

    monkeypatch.setattr("sprite_mask.validation.subprocess.run", fake_run)

    with pytest.raises(ValueError, match="do not match the sample_id: s2"):
        validate_alignment_sample_headers([Sample("s1", "popA", alignment)])


def test_validate_alignment_sample_headers_ignores_headers_without_read_group_samples(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    alignment = tmp_path / "s1.bam"
    alignment.write_bytes(b"bam")

    def fake_run(
        command: list[str],
        *,
        check: bool,
        capture_output: bool,
        text: bool,
    ) -> object:
        return subprocess.CompletedProcess(command, 0, stdout="@HD\tVN:1.6\n", stderr="")

    monkeypatch.setattr("sprite_mask.validation.subprocess.run", fake_run)

    validate_alignment_sample_headers([Sample("s1", "popA", alignment)])


def test_validate_alignment_sample_headers_reports_samtools_errors(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    alignment = tmp_path / "s1.bam"
    alignment.write_bytes(b"bam")

    def fake_run(
        command: list[str],
        *,
        check: bool,
        capture_output: bool,
        text: bool,
    ) -> object:
        return subprocess.CompletedProcess(command, 1, stdout="", stderr="not a BAM\n")

    monkeypatch.setattr("sprite_mask.validation.subprocess.run", fake_run)

    with pytest.raises(RuntimeError, match="could not read alignment header[\\s\\S]*not a BAM"):
        validate_alignment_sample_headers([Sample("s1", "popA", alignment)])


def test_require_executables_reports_missing_tools(monkeypatch: pytest.MonkeyPatch) -> None:
    def fake_which(name: str) -> str | None:
        if name == "present":
            return "/usr/bin/present"
        return None

    monkeypatch.setattr("sprite_mask.validation.shutil.which", fake_which)

    with pytest.raises(RuntimeError, match="missing1, missing2"):
        require_executables(["present", "missing1", "missing2"])


def test_ensure_parent_dirs_creates_all_parent_directories(tmp_path: Path) -> None:
    first = tmp_path / "a" / "b" / "out.bed"
    second = tmp_path / "c" / "out.bed"

    ensure_parent_dirs([first, second])

    assert first.parent.is_dir()
    assert second.parent.is_dir()


def test_refuse_existing_outputs_rejects_existing_without_force(tmp_path: Path) -> None:
    existing = tmp_path / "sprite.bed.gz"
    missing = tmp_path / "sprite.bed.gz.tbi"
    existing.write_text("old")

    with pytest.raises(FileExistsError, match=str(existing)):
        refuse_existing_outputs([existing, missing], force=False)


def test_refuse_existing_outputs_allows_existing_with_force(tmp_path: Path) -> None:
    existing = tmp_path / "sprite.bed.gz"
    existing.write_text("old")

    refuse_existing_outputs([existing], force=True)
