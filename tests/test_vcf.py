from __future__ import annotations

import gzip
import json
from pathlib import Path

import pytest

from sprite_mask.models import Sample
from sprite_mask.vcf import build_population_counts_from_all_sites_vcf, validate_vcf_sample_names


def test_build_population_counts_from_all_sites_vcf_merges_site_counts(
    tmp_path: Path,
) -> None:
    samples = [
        Sample("s1", "popA"),
        Sample("s2", "popA"),
        Sample("s3", "popB"),
    ]
    vcf = tmp_path / "all_sites.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#source\tunit-test\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\n"
        "chr1\t1\t.\tA\t.\t.\tLowQual\t.\tGT:DP\t0/0:10\t0/0:4\t0/0:8\n"
        "\n"
        "chr1\t2\t.\tC\t.\t.\t.\t.\tGT:DP\t0/0:9\t0/0:4\t0/0:1\n"
        "chr1\t2\t.\tC\tT\t.\tPASS\t.\tGT:DP\t0/0:1\t0/1:6\t0/1:10\n"
        "chr1\t3\t.\tG\t.\t.\t.\t.\tGT:DP\t0/0:0\t0/0:0\t0/0:0\n"
        "chr1\t4\t.\tT\t.\t.\t.\t.\tGT:AD:DP\t0/0:6,0:6\t0/0:5,0:5\t0/0:8,0:8\n"
        "chr1\t5\t.\tA\t.\t.\t.\t.\tGT:DP\t0/0:10\t0/0:10\t0/0:10\n"
    )
    targets = tmp_path / "targets.bed"
    targets.write_text("chr1\t0\t4\n")
    out = tmp_path / "population_counts.bed"

    build_population_counts_from_all_sites_vcf(
        samples,
        vcf,
        out,
        threshold=5,
        targets_bed=targets,
        metadata={"threshold": 5, "all_sites_vcf": str(vcf), "popfile": "popfile.tsv"},
    )

    lines = out.read_text().splitlines()
    metadata = json.loads(lines[0].split("\t", maxsplit=1)[1])
    assert metadata["format"] == "population_count_quantized_bed"
    assert metadata["source_mode"] == "all_sites_vcf"
    assert metadata["columns"] == ["chrom", "start", "end", "popA", "popB"]
    assert metadata["population_sample_counts"] == {"popA": 2, "popB": 1}
    assert lines[1:] == [
        "#chrom\tstart\tend\tpopA\tpopB",
        "chr1\t0\t1\t1\t1",
        "chr1\t1\t2\t2\t1",
        "chr1\t3\t4\t2\t1",
    ]


def test_build_population_counts_from_all_sites_vcf_requires_popfile_samples(
    tmp_path: Path,
) -> None:
    vcf = tmp_path / "all_sites.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n"
        "chr1\t1\t.\tA\t.\t.\t.\t.\tGT:DP\t0/0:10\n"
    )

    with pytest.raises(ValueError, match="popfile sample.*absent from VCF"):
        build_population_counts_from_all_sites_vcf(
            [Sample("s2", "popA")],
            vcf,
            tmp_path / "population_counts.bed",
            threshold=5,
        )


def test_build_population_counts_from_all_sites_vcf_ignores_vcf_samples_absent_from_popfile(
    tmp_path: Path,
) -> None:
    vcf = tmp_path / "all_sites.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n"
        "chr1\t1\t.\tA\t.\t.\t.\t.\tGT:DP\t0/0:0\t0/0:7\n"
    )

    out = tmp_path / "population_counts.bed"
    build_population_counts_from_all_sites_vcf(
        [Sample("s1", "popA")],
        vcf,
        out,
        threshold=5,
    )

    assert out.read_text().splitlines()[1:] == ["#chrom\tstart\tend\tpopA"]


def test_validate_vcf_sample_names_warns_for_vcf_samples_absent_from_popfile(
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    vcf = tmp_path / "all_sites.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n"
    )

    validate_vcf_sample_names([Sample("s1", "popA")], vcf)

    assert "VCF sample(s) absent from popfile" in caplog.text
    assert "these samples will be ignored: s2" in caplog.text


def test_validate_vcf_sample_names_rejects_popfile_samples_absent_from_vcf(
    tmp_path: Path,
) -> None:
    vcf = tmp_path / "all_sites.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n"
    )

    with pytest.raises(ValueError, match="popfile sample.*absent from VCF: s2"):
        validate_vcf_sample_names([Sample("s1", "popA"), Sample("s2", "popB")], vcf)


def test_build_population_counts_from_gzipped_vcf_with_custom_depth_field(
    tmp_path: Path,
) -> None:
    samples = [Sample("s1", "popA")]
    vcf = tmp_path / "all_sites.vcf.gz"
    with gzip.open(vcf, "wt") as handle:
        handle.write(
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n"
            "chr1\t1\t.\tA\t.\t.\t.\t.\tGT:MIN_DP\t0/0:7\n"
        )
    out = tmp_path / "population_counts.bed"

    build_population_counts_from_all_sites_vcf(
        samples,
        vcf,
        out,
        threshold=5,
        depth_field="MIN_DP",
    )

    assert out.read_text().splitlines()[1:] == [
        "#chrom\tstart\tend\tpopA",
        "chr1\t0\t1\t1",
    ]


def test_build_population_counts_from_all_sites_vcf_includes_non_snps_by_default(
    tmp_path: Path,
) -> None:
    samples = [Sample("s1", "popA")]
    vcf = tmp_path / "all_sites.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n"
        "chr1\t1\t.\tA\tC\t.\t.\t.\tGT:DP\t0/1:7\n"
        "chr1\t2\t.\tA\tAC\t.\t.\t.\tGT:DP\t0/1:7\n"
        "chr1\t3\t.\tAC\tA\t.\t.\t.\tGT:DP\t0/1:7\n"
        "chr1\t4\t.\tA\t.\t.\t.\t.\tGT:DP\t0/0:7\n"
        "chr1\t5\t.\tA\t<DEL>\t.\t.\t.\tGT:DP\t0/1:7\n"
        "chr1\t6\t.\tAT\tGC\t.\t.\t.\tGT:DP\t0/1:7\n"
        "chr1\t7\t.\tA\tC,G\t.\t.\t.\tGT:DP\t1/2:7\n"
        "chr1\t8\t.\tA\tC,AC\t.\t.\t.\tGT:DP\t1/2:7\n"
    )
    out = tmp_path / "population_counts.bed"

    build_population_counts_from_all_sites_vcf(samples, vcf, out, threshold=5)

    assert out.read_text().splitlines()[1:] == [
        "#chrom\tstart\tend\tpopA",
        "chr1\t0\t8\t1",
    ]


def test_build_population_counts_from_all_sites_vcf_snps_only_keeps_snps_and_invariant_sites(
    tmp_path: Path,
) -> None:
    samples = [Sample("s1", "popA")]
    vcf = tmp_path / "all_sites.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n"
        "chr1\t1\t.\tA\tC\t.\t.\t.\tGT:DP\t0/1:7\n"
        "chr1\t2\t.\tA\tAC\t.\t.\t.\tGT:DP\t0/1:7\n"
        "chr1\t3\t.\tAC\tA\t.\t.\t.\tGT:DP\t0/1:7\n"
        "chr1\t4\t.\tA\t.\t.\t.\t.\tGT:DP\t0/0:7\n"
        "chr1\t5\t.\tA\t<DEL>\t.\t.\t.\tGT:DP\t0/1:7\n"
        "chr1\t6\t.\tAT\tGC\t.\t.\t.\tGT:DP\t0/1:7\n"
        "chr1\t7\t.\tA\tC,G\t.\t.\t.\tGT:DP\t1/2:7\n"
        "chr1\t8\t.\tA\tC,AC\t.\t.\t.\tGT:DP\t1/2:7\n"
    )
    out = tmp_path / "population_counts.bed"

    build_population_counts_from_all_sites_vcf(samples, vcf, out, threshold=5, snps_only=True)

    assert out.read_text().splitlines()[1:] == [
        "#chrom\tstart\tend\tpopA",
        "chr1\t0\t1\t1",
        "chr1\t3\t4\t1",
        "chr1\t6\t7\t1",
    ]


def test_build_population_counts_from_all_sites_vcf_snps_only_uses_snp_duplicate_records(
    tmp_path: Path,
) -> None:
    samples = [Sample("s1", "popA")]
    vcf = tmp_path / "all_sites.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n"
        "chr1\t1\t.\tA\tAC\t.\t.\t.\tGT:DP\t0/1:7\n"
        "chr1\t1\t.\tA\tC\t.\t.\t.\tGT:DP\t0/1:7\n"
        "chr1\t2\t.\tA\t.\t.\t.\t.\tGT:DP\t0/0:7\n"
    )
    out = tmp_path / "population_counts.bed"

    build_population_counts_from_all_sites_vcf(samples, vcf, out, threshold=5, snps_only=True)

    assert out.read_text().splitlines()[1:] == [
        "#chrom\tstart\tend\tpopA",
        "chr1\t0\t2\t1",
    ]


def test_build_population_counts_from_all_sites_vcf_respects_target_boundaries(
    tmp_path: Path,
) -> None:
    samples = [Sample("s1", "popA")]
    vcf = tmp_path / "all_sites.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n"
        "chr1\t1\t.\tA\t.\t.\t.\t.\tGT:DP\t0/0:7\n"
        "chr1\t2\t.\tC\t.\t.\t.\t.\tGT:DP\t0/0:7\n"
        "chr1\t3\t.\tG\t.\t.\t.\t.\tGT:DP\t0/0:7\n"
        "chr1\t4\t.\tT\t.\t.\t.\t.\tGT:DP\t0/0:7\n"
        "chr1\t5\t.\tA\t.\t.\t.\t.\tGT:DP\t0/0:7\n"
        "chr2\t1\t.\tA\t.\t.\t.\t.\tGT:DP\t0/0:7\n"
    )
    targets = tmp_path / "targets.bed"
    targets.write_text("chr1\t1\t2\nchr1\t2\t3\nchr1\t4\t5\n")
    out = tmp_path / "population_counts.bed"

    build_population_counts_from_all_sites_vcf(
        samples,
        vcf,
        out,
        threshold=5,
        targets_bed=targets,
    )

    assert out.read_text().splitlines()[1:] == [
        "#chrom\tstart\tend\tpopA",
        "chr1\t1\t3\t1",
        "chr1\t4\t5\t1",
    ]


def test_build_population_counts_from_all_sites_vcf_rejects_noncontiguous_duplicates(
    tmp_path: Path,
) -> None:
    samples = [Sample("s1", "popA")]
    vcf = tmp_path / "all_sites.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n"
        "chr1\t1\t.\tA\t.\t.\t.\t.\tGT:DP\t0/0:0\n"
        "chr1\t2\t.\tC\t.\t.\t.\t.\tGT:DP\t0/0:7\n"
        "chr1\t1\t.\tA\tT\t.\t.\t.\tGT:DP\t0/1:7\n"
    )

    with pytest.raises(ValueError, match="duplicate CHROM:POS records are contiguous"):
        build_population_counts_from_all_sites_vcf(
            samples,
            vcf,
            tmp_path / "population_counts.bed",
            threshold=5,
        )


def test_build_population_counts_from_all_sites_vcf_rejects_reopened_chromosomes(
    tmp_path: Path,
) -> None:
    samples = [Sample("s1", "popA")]
    vcf = tmp_path / "all_sites.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n"
        "chr1\t1\t.\tA\t.\t.\t.\t.\tGT:DP\t0/0:7\n"
        "chr2\t1\t.\tC\t.\t.\t.\t.\tGT:DP\t0/0:7\n"
        "chr1\t2\t.\tG\t.\t.\t.\t.\tGT:DP\t0/0:7\n"
    )

    with pytest.raises(ValueError, match="grouped by chromosome"):
        build_population_counts_from_all_sites_vcf(
            samples,
            vcf,
            tmp_path / "population_counts.bed",
            threshold=5,
        )


def test_build_population_counts_from_all_sites_vcf_omits_missing_and_zero_depths(
    tmp_path: Path,
) -> None:
    samples = [Sample("s1", "popA")]
    vcf = tmp_path / "all_sites.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n"
        "chr1\t1\t.\tA\t.\t.\t.\t.\tGT:DP\t.\n"
        "chr1\t2\t.\tC\t.\t.\t.\t.\tGT:DP\t0/0:.\n"
        "chr1\t3\t.\tG\t.\t.\t.\t.\tGT:DP\t0/0\n"
        "chr1\t4\t.\tT\t.\t.\t.\t.\tGT:DP\t0/0:0\n"
    )
    out = tmp_path / "population_counts.bed"

    build_population_counts_from_all_sites_vcf(samples, vcf, out, threshold=5)

    assert out.read_text().splitlines()[1:] == ["#chrom\tstart\tend\tpopA"]


def test_build_population_counts_from_all_sites_vcf_with_no_records_writes_only_header(
    tmp_path: Path,
) -> None:
    samples = [Sample("s1", "popA")]
    vcf = tmp_path / "all_sites.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n"
    )
    out = tmp_path / "population_counts.bed"

    build_population_counts_from_all_sites_vcf(samples, vcf, out, threshold=5)

    assert out.read_text().splitlines()[1:] == ["#chrom\tstart\tend\tpopA"]


@pytest.mark.parametrize(
    ("contents", "message"),
    [
        (
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n",
            "does not contain VCF sample columns",
        ),
        (
            "chr1\t1\t.\tA\t.\t.\t.\t.\tGT:DP\t0/0:7\n",
            "appears before the #CHROM header",
        ),
        (
            "##fileformat=VCFv4.2\n",
            "does not contain a #CHROM header",
        ),
        (
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts1\n",
            "duplicate VCF sample",
        ),
    ],
)
def test_build_population_counts_from_all_sites_vcf_rejects_malformed_headers(
    tmp_path: Path,
    contents: str,
    message: str,
) -> None:
    vcf = tmp_path / "all_sites.vcf"
    vcf.write_text(contents)

    with pytest.raises(ValueError, match=message):
        build_population_counts_from_all_sites_vcf(
            [Sample("s1", "popA")],
            vcf,
            tmp_path / "population_counts.bed",
            threshold=5,
        )


@pytest.mark.parametrize(
    ("header_samples", "samples", "row", "message"),
    [
        (
            ["s1"],
            [Sample("s1", "popA")],
            "chr1\t1\t.\tA\t.\t.\t.\t.\n",
            "must have VCF fixed fields",
        ),
        (
            ["s1"],
            [Sample("s1", "popA")],
            "chr1\tbad\t.\tA\t.\t.\t.\t.\tGT:DP\t0/0:7\n",
            "non-integer POS",
        ),
        (
            ["s1"],
            [Sample("s1", "popA")],
            "chr1\t0\t.\tA\t.\t.\t.\t.\tGT:DP\t0/0:7\n",
            "POS < 1",
        ),
        (
            ["s1"],
            [Sample("s1", "popA")],
            "chr1\t1\t.\tA\t.\t.\t.\t.\tGT:AD\t0/0:7\n",
            "FORMAT does not contain",
        ),
        (
            ["s1", "s2"],
            [Sample("s1", "popA"), Sample("s2", "popA")],
            "chr1\t1\t.\tA\t.\t.\t.\t.\tGT:DP\t0/0:7\n",
            "fewer sample columns",
        ),
        (
            ["s1"],
            [Sample("s1", "popA")],
            "chr1\t1\t.\tA\t.\t.\t.\t.\tGT:DP\t0/0:nope\n",
            "non-integer sample DP",
        ),
    ],
)
def test_build_population_counts_from_all_sites_vcf_rejects_malformed_records(
    tmp_path: Path,
    header_samples: list[str],
    samples: list[Sample],
    row: str,
    message: str,
) -> None:
    vcf = tmp_path / "all_sites.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        + "\t".join(header_samples)
        + "\n"
        + row
    )

    with pytest.raises(ValueError, match=message):
        build_population_counts_from_all_sites_vcf(
            samples,
            vcf,
            tmp_path / "population_counts.bed",
            threshold=5,
        )
