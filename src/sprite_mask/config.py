from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class AlignmentRunConfig:
    samples_path: Path
    min_dp: int
    out_dir: Path
    work_dir: Path | None = None
    threads: int = 1
    jobs: int = 1
    mask_bed: Path | None = None
    min_mapq: int | None = None
    max_dp: int | None = None
    exclude_flag: int | None = None
    reference: Path | None = None
    strict_depth: bool = False
    keep_work: bool = False
    force: bool = False
    dry_run: bool = False
    output_prefix: str = "sprite"

    @property
    def resolved_work_dir(self) -> Path:
        return self.work_dir if self.work_dir is not None else self.out_dir / "work"


@dataclass(frozen=True)
class VcfRunConfig:
    all_sites_vcf: Path
    popfile_path: Path
    min_dp: int
    out_dir: Path
    work_dir: Path | None = None
    mask_bed: Path | None = None
    keep_work: bool = False
    force: bool = False
    dry_run: bool = False
    output_prefix: str = "sprite"
    snps_only: bool = False

    @property
    def resolved_work_dir(self) -> Path:
        return self.work_dir if self.work_dir is not None else self.out_dir / "work"


RunConfig = AlignmentRunConfig | VcfRunConfig
