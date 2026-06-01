from __future__ import annotations

import subprocess
from collections.abc import Sequence
from pathlib import Path

from sprite_mask.bedio import iter_bed3
from sprite_mask.commands import run_pipeline


def sort_and_merge_bed(in_bed: Path, out_bed: Path) -> Path:
    run_pipeline(
        [
            ["bedtools", "sort", "-i", str(in_bed)],
            ["bedtools", "merge", "-i", "-"],
        ],
        out_bed,
    )
    return out_bed


def intersect_sort_merge(a_bed: Path, b_bed: Path, out_bed: Path) -> Path:
    run_pipeline(
        [
            ["bedtools", "intersect", "-a", str(a_bed), "-b", str(b_bed)],
            ["bedtools", "sort", "-i", "-"],
            ["bedtools", "merge", "-i", "-"],
        ],
        out_bed,
    )
    return out_bed


def subtract_sort_merge(a_bed: Path, b_bed: Path, out_bed: Path) -> Path:
    run_pipeline(
        [
            ["bedtools", "subtract", "-a", str(a_bed), "-b", str(b_bed)],
            ["bedtools", "sort", "-i", "-"],
            ["bedtools", "merge", "-i", "-"],
        ],
        out_bed,
    )
    return out_bed


def run_multiinter(pass_beds: Sequence[Path], names: Sequence[str], out_tsv: Path) -> Path:
    command = build_multiinter_command(pass_beds, names)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with out_tsv.open("w") as out:
        subprocess.run(command, check=True, stdout=out, text=True)
    return out_tsv


def write_single_input_multiinter(pass_bed: Path, name: str, out_tsv: Path) -> Path:
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with out_tsv.open("w") as out:
        out.write(f"chrom\tstart\tend\tnum\tlist\t{name}\n")
        for chrom, start, end in iter_bed3(pass_bed):
            out.write(f"{chrom}\t{start}\t{end}\t1\t{name}\t1\n")
    return out_tsv


def build_multiinter_command(pass_beds: Sequence[Path], names: Sequence[str]) -> list[str]:
    if not pass_beds:
        raise ValueError("at least one pass BED is required")
    if len(pass_beds) != len(names):
        raise ValueError("pass BED count must match sample name count")
    return [
        "bedtools",
        "multiinter",
        "-header",
        "-names",
        *names,
        "-i",
        *[str(path) for path in pass_beds],
    ]
