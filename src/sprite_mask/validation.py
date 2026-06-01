from __future__ import annotations

import shutil
import subprocess
from collections.abc import Iterable
from pathlib import Path

from sprite_mask.models import Sample


def validate_threshold(
    threshold: int,
    *,
    targets_bed: Path | None = None,
    all_sites_vcf: Path | None = None,
) -> None:
    if threshold < 0:
        raise ValueError("--threshold must be a non-negative integer")
    if threshold == 0 and targets_bed is None and all_sites_vcf is None:
        raise ValueError("--min-dp 0 requires --mask because no genome.txt is required")


def validate_threads(threads: int) -> None:
    if threads < 1:
        raise ValueError("--threads must be at least 1")


def validate_jobs(jobs: int) -> None:
    if jobs < 1:
        raise ValueError("--jobs must be at least 1")


def validate_vcf_inputs(all_sites_vcf: Path, popfile_path: Path) -> None:
    if not all_sites_vcf.exists():
        raise ValueError(f"--all-sites-vcf does not exist: {all_sites_vcf}")
    if all_sites_vcf.is_dir():
        raise ValueError(f"--all-sites-vcf is a directory: {all_sites_vcf}")
    if not popfile_path.exists():
        raise ValueError(f"--popfile does not exist: {popfile_path}")
    if popfile_path.is_dir():
        raise ValueError(f"--popfile is a directory: {popfile_path}")


def validate_variants_vcf_input(variants_vcf: Path) -> None:
    if not variants_vcf.exists():
        raise ValueError(f"--variants-vcf does not exist: {variants_vcf}")
    if variants_vcf.is_dir():
        raise ValueError(f"--variants-vcf is a directory: {variants_vcf}")


def validate_alignment_sample_headers(samples: list[Sample]) -> None:
    for sample in samples:
        if sample.alignment is None:
            raise ValueError(f"alignment for sample {sample.sample_id!r} is missing")

        header = _read_alignment_header(sample.alignment, sample.sample_id)
        header_sample_names = _read_group_sample_names(header)
        unexpected_sample_names = sorted(
            sample_name for sample_name in header_sample_names if sample_name != sample.sample_id
        )
        if unexpected_sample_names:
            raise ValueError(
                f"alignment for sample {sample.sample_id!r} contains @RG SM value(s) "
                f"that do not match the sample_id: {', '.join(unexpected_sample_names)} "
                f"({sample.alignment})"
            )


def _read_alignment_header(alignment: Path, sample_id: str) -> str:
    completed = subprocess.run(
        ["samtools", "view", "-H", str(alignment)],
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        message = completed.stderr.strip() or completed.stdout.strip() or "unknown samtools error"
        raise RuntimeError(
            f"could not read alignment header for sample {sample_id!r}: {alignment}\n{message}"
        )
    return completed.stdout


def _read_group_sample_names(header: str) -> set[str]:
    sample_names: set[str] = set()
    for line in header.splitlines():
        if not line.startswith("@RG\t"):
            continue
        for field in line.split("\t")[1:]:
            if field.startswith("SM:") and len(field) > 3:
                sample_names.add(field[3:])
    return sample_names


def require_executables(names: Iterable[str]) -> None:
    missing = [name for name in names if shutil.which(name) is None]
    if missing:
        raise RuntimeError(f"required executable(s) not found on PATH: {', '.join(missing)}")


def ensure_parent_dirs(paths: Iterable[Path]) -> None:
    for path in paths:
        path.parent.mkdir(parents=True, exist_ok=True)


def refuse_existing_outputs(paths: Iterable[Path], *, force: bool) -> None:
    existing = [path for path in paths if path.exists()]
    if existing and not force:
        joined = "\n".join(str(path) for path in existing)
        raise FileExistsError(f"output file(s) already exist; pass --force to overwrite:\n{joined}")
