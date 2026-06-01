from __future__ import annotations

import gzip
import json
import shutil
import subprocess
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
TEST_THRESHOLD = 30

CountBedRow = tuple[str, int, int, tuple[int, ...]]
TargetInterval = tuple[str, int, int]


@dataclass(frozen=True)
class FixtureCase:
    name: str
    fixture: Path
    target_intervals: tuple[TargetInterval, ...]
    expected_sample_counts: dict[str, int]
    expected_sample_summary: dict[int, int]
    expected_population_summary: dict[str, dict[int, int]]
    expected_population_counts: dict[str, int]
    min_count_rows: int

    @property
    def samples_tsv(self) -> Path:
        return self.fixture / "samples.tsv"

    @property
    def targets_bed(self) -> Path:
        return self.fixture / "targets.bed"

    @property
    def populations(self) -> list[str]:
        return list(self.expected_population_counts)

    @property
    def target_sites(self) -> int:
        return sum(end - start for _chrom, start, end in self.target_intervals)


SMOKE_FIXTURE = FixtureCase(
    name="1000g_5sample_chr20_smoke",
    fixture=REPO_ROOT / "tests" / "test_data" / "1000g_5sample_chr20_smoke",
    target_intervals=(("chr20", 9_999_999, 10_009_999),),
    expected_sample_counts={
        "HG00096": 0,
        "HG00097": 0,
        "HG00099": 0,
        "NA18486": 0,
        "NA18489": 0,
    },
    expected_sample_summary={},
    expected_population_summary={},
    expected_population_counts={"GBR": 3, "YRI": 2},
    min_count_rows=1,
)

INTEGRATION_FIXTURE = FixtureCase(
    name="1000g_20sample_highcov_4chrom_subset",
    fixture=REPO_ROOT / "tests" / "test_data" / "1000g_20sample_highcov_4chrom_subset",
    target_intervals=(
        ("chr2", 9_999_999, 10_049_999),
        ("chr7", 49_999_999, 50_049_999),
        ("chr12", 24_999_999, 25_049_999),
        ("chr20", 9_999_999, 10_049_999),
    ),
    expected_sample_counts={
        "HG00096": 27_521,
        "HG00097": 20_585,
        "HG00099": 58_048,
        "HG00100": 14_469,
        "HG00101": 28_605,
        "HG00102": 18_908,
        "HG00103": 22_460,
        "HG00105": 8_039,
        "HG00106": 47_679,
        "HG00107": 22_583,
        "NA18486": 21_089,
        "NA18488": 10_702,
        "NA18489": 190_337,
        "NA18498": 8_661,
        "NA18499": 24_514,
        "NA18501": 12_409,
        "NA18502": 9_021,
        "NA18504": 12_821,
        "NA18505": 50_961,
        "NA18507": 24_492,
    },
    expected_sample_summary={
        1: 36_250,
        2: 49_967,
        3: 42_286,
        4: 28_271,
        5: 17_221,
        6: 9_096,
        7: 4_991,
        8: 3_030,
        9: 1_822,
        10: 1_134,
        11: 693,
        12: 439,
        13: 291,
        14: 171,
        15: 96,
        16: 76,
        17: 178,
        18: 129,
        19: 130,
        20: 32,
    },
    expected_population_summary={
        "GBR": {
            0: 57_972,
            1: 64_395,
            2: 40_867,
            3: 19_713,
            4: 7_959,
            5: 2_947,
            6: 1_177,
            7: 617,
            8: 299,
            9: 285,
            10: 72,
        },
        "YRI": {
            0: 2_579,
            1: 90_131,
            2: 62_499,
            3: 25_836,
            4: 9_429,
            5: 3_091,
            6: 1_252,
            7: 789,
            8: 232,
            9: 342,
            10: 123,
        },
    },
    expected_population_counts={"GBR": 10, "YRI": 10},
    min_count_rows=30_000,
)


FAST_MODE_INTEGRATION_FIXTURE = FixtureCase(
    name="1000g_20sample_highcov_4chrom_subset",
    fixture=REPO_ROOT / "tests" / "test_data" / "1000g_20sample_highcov_4chrom_subset",
    target_intervals=INTEGRATION_FIXTURE.target_intervals,
    expected_sample_counts={
        "HG00096": 29_770,
        "HG00097": 22_129,
        "HG00099": 61_681,
        "HG00100": 16_207,
        "HG00101": 30_339,
        "HG00102": 20_120,
        "HG00103": 24_599,
        "HG00105": 9_139,
        "HG00106": 50_536,
        "HG00107": 24_792,
        "NA18486": 23_563,
        "NA18488": 11_644,
        "NA18489": 191_271,
        "NA18498": 10_303,
        "NA18499": 26_477,
        "NA18501": 13_995,
        "NA18502": 10_405,
        "NA18504": 14_134,
        "NA18505": 54_741,
        "NA18507": 26_329,
    },
    expected_sample_summary={
        1: 31_321,
        2: 47_522,
        3: 42_668,
        4: 30_351,
        5: 19_167,
        6: 10_652,
        7: 5_799,
        8: 3_505,
        9: 2_040,
        10: 1_297,
        11: 839,
        12: 509,
        13: 394,
        14: 209,
        15: 108,
        16: 93,
        17: 116,
        18: 180,
        19: 174,
        20: 84,
    },
    expected_population_summary={
        "GBR": {
            0: 52_576,
            1: 64_140,
            2: 43_232,
            3: 21_653,
            4: 9_160,
            5: 3_460,
            6: 1_356,
            7: 654,
            8: 310,
            9: 255,
            10: 232,
        },
        "YRI": {
            0: 2_326,
            1: 83_791,
            2: 64_760,
            3: 28_526,
            4: 10_859,
            5: 3_483,
            6: 1_506,
            7: 859,
            8: 366,
            9: 375,
            10: 177,
        },
    },
    expected_population_counts=INTEGRATION_FIXTURE.expected_population_counts,
    min_count_rows=INTEGRATION_FIXTURE.min_count_rows,
)


def test_cli_workflow_multichromosome_fixture_outputs_expected_counts(tmp_path: Path) -> None:
    fixture_case = INTEGRATION_FIXTURE
    require_fixture_and_tools(fixture_case)
    out_dir, work_dir = run_sprite_mask(tmp_path, fixture_case, keep_work=True)

    assert_population_count_output_matches_fixture(out_dir, fixture_case)
    assert_expected_work_files_exist(work_dir, fixture_case)


def test_cli_workflow_fast_mode_preserves_previous_fixture_counts(tmp_path: Path) -> None:
    fixture_case = FAST_MODE_INTEGRATION_FIXTURE
    require_fixture_and_tools(fixture_case)
    out_dir, work_dir = run_sprite_mask(
        tmp_path,
        fixture_case,
        keep_work=True,
        fast_mode=True,
    )

    assert_population_count_output_matches_fixture(out_dir, fixture_case)
    assert_expected_work_files_exist(work_dir, fixture_case)


def assert_population_count_output_matches_fixture(
    out_dir: Path,
    fixture_case: FixtureCase,
) -> None:
    population_count_bed_gz = out_dir / "sprite.bed.gz"
    population_count_bed_index = Path(f"{population_count_bed_gz}.tbi")

    assert sorted(path.name for path in out_dir.iterdir()) == [
        population_count_bed_gz.name,
        population_count_bed_index.name,
    ]
    assert population_count_bed_gz.stat().st_size > 0
    assert population_count_bed_index.stat().st_size > 0

    population_rows = read_count_bed_gz(
        population_count_bed_gz,
        expected_header=fixture_case.populations,
    )
    population_metadata = read_count_bed_header(population_count_bed_gz)
    assert population_metadata["columns"] == ["chrom", "start", "end", *fixture_case.populations]
    assert population_metadata["populations"] == fixture_case.populations
    assert query_tabix_rows(population_count_bed_gz, fixture_case.target_intervals[0]) > 0

    assert len(population_rows) > fixture_case.min_count_rows
    assert summarize_population_bed(population_rows, fixture_case.populations) == (
        fixture_case.expected_population_summary
    )
    assert summarize_population_totals(population_rows) == fixture_case.expected_sample_summary

    expected_nonzero_sites = sum(fixture_case.expected_sample_summary.values())
    assert_count_bed_shape(
        population_rows,
        target_intervals=fixture_case.target_intervals,
        value_count=2,
        max_values=list(fixture_case.expected_population_counts.values()),
        expected_covered_sites=expected_nonzero_sites,
        expected_target_sites=fixture_case.target_sites,
    )


def test_cli_workflow_smoke_fixture_writes_only_indexed_population_bed(tmp_path: Path) -> None:
    fixture_case = SMOKE_FIXTURE
    require_fixture_and_tools(fixture_case)
    out_dir, work_dir = run_sprite_mask(tmp_path, fixture_case, keep_work=False, threads=1)

    population_count_bed_gz = out_dir / "sprite.bed.gz"
    output_names = sorted(path.name for path in out_dir.iterdir())
    assert output_names == [
        population_count_bed_gz.name,
        f"{population_count_bed_gz.name}.tbi",
    ]

    metadata = read_count_bed_header(population_count_bed_gz)
    assert metadata["columns"] == ["chrom", "start", "end", *fixture_case.populations]
    assert metadata["populations"] == fixture_case.populations
    assert metadata["population_sample_counts"] == fixture_case.expected_population_counts
    assert query_tabix_rows(population_count_bed_gz, fixture_case.target_intervals[0]) > 0
    assert not work_dir.exists()


def require_fixture_and_tools(fixture_case: FixtureCase) -> None:
    if not fixture_case.samples_tsv.exists() or not fixture_case.targets_bed.exists():
        pytest.skip(f"1000 Genomes BAM fixture is not available: {fixture_case.name}")
    missing = [
        tool
        for tool in ("samtools", "mosdepth", "bedtools", "bgzip", "tabix")
        if shutil.which(tool) is None
    ]
    if missing:
        pytest.skip(f"external workflow tool(s) unavailable: {', '.join(missing)}")


def run_sprite_mask(
    tmp_path: Path,
    fixture_case: FixtureCase,
    *,
    keep_work: bool,
    threads: int = 2,
    jobs: int = 2,
    fast_mode: bool = False,
) -> tuple[Path, Path]:
    out_dir = tmp_path / "results"
    work_dir = tmp_path / "work"
    command = [
        sys.executable,
        "-m",
        "sprite_mask.cli",
        "from-alignments",
        "--samples",
        str(fixture_case.samples_tsv),
        "--min-dp",
        str(TEST_THRESHOLD),
        "--mask",
        str(fixture_case.targets_bed),
        "--out",
        str(out_dir),
        "--work",
        str(work_dir),
        "--threads",
        str(threads),
        "--jobs",
        str(jobs),
    ]
    if fast_mode:
        command.append("--fast-mode")
    if keep_work:
        command.append("--keep-work")

    completed = subprocess.run(
        command,
        cwd=REPO_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    assert completed.stdout == ""
    assert "[sprite] Analysis complete: wrote" in completed.stderr
    assert "sprite.bed.gz" in completed.stderr
    return out_dir, work_dir


def read_count_bed_gz(
    path: Path,
    *,
    expected_header: list[str],
) -> list[CountBedRow]:
    with gzip.open(path, "rt") as handle:
        header = _read_count_bed_column_header(handle)
        assert header == ["chrom", "start", "end", *expected_header]
        rows = []
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            rows.append((fields[0], int(fields[1]), int(fields[2]), tuple(map(int, fields[3:]))))
    return rows


def read_count_bed_header(path: Path) -> dict[str, object]:
    with gzip.open(path, "rt") as handle:
        metadata: dict[str, object] | None = None
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields[0] == "#sprite_mask_metadata":
                metadata = json.loads(fields[1])
                continue
            if fields[0] == "#chrom":
                assert metadata is not None
                assert [fields[0].lstrip("#"), *fields[1:]] == metadata["columns"]
                return metadata
    raise AssertionError(f"{path} does not contain a quantized BED header")


def _read_count_bed_column_header(handle: object) -> list[str]:
    for line in handle:
        fields = line.rstrip("\n").split("\t")
        if fields[0] == "#chrom":
            return [fields[0].lstrip("#"), *fields[1:]]
        if fields[0] == "chrom":
            return fields
    raise AssertionError("count BED does not contain a column header")


def query_tabix_rows(path: Path, target_interval: TargetInterval) -> int:
    chrom, start, end = target_interval
    completed = subprocess.run(
        ["tabix", str(path), f"{chrom}:{start + 1}-{end}"],
        check=True,
        capture_output=True,
        text=True,
    )
    return len(completed.stdout.splitlines())


def assert_count_bed_shape(
    rows: list[CountBedRow],
    *,
    target_intervals: tuple[TargetInterval, ...],
    value_count: int,
    max_values: list[int],
    expected_covered_sites: int,
    expected_target_sites: int,
) -> None:
    previous: CountBedRow | None = None
    covered_sites = 0
    for chrom, start, end, values in rows:
        assert any(
            target_chrom == chrom and target_start <= start < end <= target_end
            for target_chrom, target_start, target_end in target_intervals
        )
        assert len(values) == value_count
        assert all(
            0 <= value <= max_value
            for value, max_value in zip(values, max_values, strict=True)
        )
        covered_sites += end - start

        if previous is not None:
            previous_chrom, _previous_start, previous_end, previous_values = previous
            assert (previous_chrom, previous_end) <= (chrom, start)
            if previous_chrom == chrom and previous_end == start:
                assert previous_values != values
        previous = (chrom, start, end, values)

    assert 0 < covered_sites <= expected_target_sites
    assert covered_sites == expected_covered_sites


def summarize_population_bed(
    rows: list[CountBedRow],
    populations: list[str],
) -> dict[str, dict[int, int]]:
    summary: dict[str, Counter[int]] = {population: Counter() for population in populations}
    for _chrom, start, end, values in rows:
        length = end - start
        for population, value in zip(populations, values, strict=True):
            summary[population][value] += length
    return {population: dict(counts) for population, counts in summary.items()}


def summarize_population_totals(
    rows: list[CountBedRow],
) -> dict[int, int]:
    summary: Counter[int] = Counter()
    for _chrom, start, end, values in rows:
        summary[sum(values)] += end - start
    return dict(summary)


def assert_expected_work_files_exist(work_dir: Path, fixture_case: FixtureCase) -> None:
    assert work_dir.exists()
    expected_sample_suffixes = (
        "quantized.bed.gz",
        "quantized.bed.gz.csi",
        "mosdepth.summary.txt",
        "mosdepth.global.dist.txt",
        "mosdepth.stderr.log",
        "pass.bed",
        "pass.targets.bed",
    )
    for sample_id in fixture_case.expected_sample_counts:
        for suffix in expected_sample_suffixes:
            assert (work_dir / f"{sample_id}.d{TEST_THRESHOLD}.{suffix}").exists()

    assert (work_dir / "targets.3col.bed").exists()
    assert (work_dir / "targets.3col.sorted.merged.bed").exists()
    assert (work_dir / f"cohort.d{TEST_THRESHOLD}.multiinter.tsv").exists()
    assert (work_dir / f"cohort.d{TEST_THRESHOLD}.population_count_quantized.bed").exists()
