from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from sprite_mask.summaries import summarize_population_count_bed


def test_summarize_population_count_bed(tmp_path: Path) -> None:
    bed = tmp_path / "population_counts.bed"
    bed.write_text(
        "#sprite_mask_metadata\t{}\n"
        "#chrom\tstart\tend\tpopA\tpopB\n"
        "chr1\t0\t10\t1\t0\n"
        "chr1\t10\t25\t2\t1\n"
        "chr2\t0\t5\t1\t0\n"
    )
    out = tmp_path / "summary.tsv"

    summarize_population_count_bed(bed, out)

    assert out.read_text() == (
        "population\tpassing_samples\tsites\n"
        "popA\t1\t15\n"
        "popA\t2\t15\n"
        "popB\t0\t15\n"
        "popB\t1\t15\n"
    )


def test_summarize_population_count_bed_accepts_uncommented_header(
    tmp_path: Path,
) -> None:
    bed = tmp_path / "population_counts.bed"
    bed.write_text(
        "chrom\tstart\tend\tpopA\n"
        "chr1\t0\t10\t1\n"
        "chr1\t10\t15\t0\n"
    )
    out = tmp_path / "summary.tsv"

    summarize_population_count_bed(bed, out)

    assert out.read_text() == (
        "population\tpassing_samples\tsites\n"
        "popA\t0\t5\n"
        "popA\t1\t10\n"
    )


def test_summarize_population_count_bed_accepts_gzipped_input(tmp_path: Path) -> None:
    bed = tmp_path / "population_counts.bed.gz"
    with gzip.open(bed, "wt") as handle:
        handle.write(
            "#chrom\tstart\tend\tpopA\n"
            "chr1\t0\t10\t1\n"
            "chr1\t10\t15\t0\n"
        )
    out = tmp_path / "summary.tsv"

    summarize_population_count_bed(bed, out)

    assert out.read_text() == (
        "population\tpassing_samples\tsites\n"
        "popA\t0\t5\n"
        "popA\t1\t10\n"
    )


def test_summarize_population_count_bed_rejects_empty_input(tmp_path: Path) -> None:
    bed = tmp_path / "population_counts.bed"
    bed.write_text("")

    with pytest.raises(ValueError, match="is empty"):
        summarize_population_count_bed(bed, tmp_path / "summary.tsv")


def test_summarize_population_count_bed_rejects_missing_header(tmp_path: Path) -> None:
    bed = tmp_path / "population_counts.bed"
    bed.write_text("chr1\t0\t10\t1\n")

    with pytest.raises(ValueError, match="must have a header"):
        summarize_population_count_bed(bed, tmp_path / "summary.tsv")


def test_summarize_population_count_bed_skips_blank_and_comment_lines(
    tmp_path: Path,
) -> None:
    bed = tmp_path / "population_counts.bed"
    bed.write_text(
        "\n"
        "#not_the_header\tignored\n"
        "#chrom\tstart\tend\tpopA\n"
        "\n"
        "# comment\n"
        "chrom\tstart\tend\tpopA\n"
        "chr1\t0\t3\t1\n"
    )
    out = tmp_path / "summary.tsv"

    summarize_population_count_bed(bed, out)

    assert out.read_text() == "population\tpassing_samples\tsites\npopA\t1\t3\n"
