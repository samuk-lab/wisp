# Implementation plan

## Goal

Build a single Python CLI tool that calculates depth-threshold pass intervals per BAM/CRAM with `mosdepth`, intersects those intervals across samples with `bedtools multiinter`, and writes collapsed quantized BED outputs:

- Per-sample counts: number of sites in each BAM/CRAM that meet the depth threshold.
- Cohort count BED: intervals labelled by the number of samples that meet the threshold.
- Population count BED: intervals with one count column per population.
- Summary TSVs for sample-count and population-count distributions.

The tool should use sparse output: intervals absent from the final BED are assumed to have zero passing samples. Do not require users to pass a `genome.txt` file.

## Package Shape

Use a Python package with one primary command, plus small debug/development subcommands if useful.

```text
sprite-mask/
  environment.yml
  pyproject.toml
  README.md
  src/sprite_mask/
    __init__.py
    cli.py
    config.py
    models.py
    samples.py
    workflow.py
    commands.py
    mosdepth.py
    bedtools.py
    bedio.py
    collapse.py
    summaries.py
    validation.py
  tests/
    data/
    test_samples.py
    test_bedio.py
    test_collapse.py
    test_summaries.py
    test_workflow_commands.py
```

Use `pyproject.toml` with a console entry point:

```toml
[project.scripts]
sprite-mask = "sprite_mask.cli:main"
```

## Development Environment

Development can happen in either macOS or WSL, depending on where the current Codex session is running. Use the same conda environment name, `sprite`, and the same `environment.yml` in both places.

When working under WSL, prefer keeping the checkout on the WSL/Linux filesystem, such as `~/projects/sprite-mask`, rather than under `/mnt/c`, because the workflow will create and read many intermediate files. When working under macOS, use the local macOS checkout and the macOS `sprite` conda environment.

Use the mamba/libmamba solver when possible:

```bash
mamba env create -f environment.yml
conda activate sprite
```

If the environment already exists:

```bash
mamba env update -n sprite -f environment.yml
conda activate sprite
```

If `mamba` is not installed, use conda with the libmamba solver when available:

```bash
conda env create -f environment.yml --solver=libmamba
conda activate sprite
```

After the package skeleton exists, install the project in editable mode:

```bash
python -m pip install -e ".[dev]"
```

Starter `environment.yml`:

```yaml
name: sprite
channels:
  - conda-forge
  - bioconda
dependencies:
  - python>=3.11,<3.13
  - pip
  - mosdepth
  - bedtools
  - samtools
  - pytest
  - pytest-cov
  - ruff
  - mypy
  - build
```

## Build And CI Targets

Every build should run on both Linux and macOS.

Use a CI matrix with:

- `ubuntu-latest`
- `macos-latest`

Each matrix job should:

1. Check out the repository.
2. Create or update the `sprite` conda environment from `environment.yml`, using mamba/micromamba where possible.
3. Activate `sprite`.
4. Install the package in editable mode with development extras.
5. Run linting and tests.
6. Build the Python source distribution and wheel.

Command shape inside each CI job:

```bash
python -m pip install -e ".[dev]"
ruff check .
mypy src
pytest
python -m build
```

Because this is planned as a pure Python package that orchestrates external binaries, the wheel should be platform independent. The reason to build and test on both Linux and macOS is to verify subprocess behavior, paths, conda dependency resolution, and `mosdepth`/`bedtools` availability on both platforms.

The CI environment should use the same dependency channels as local development:

```text
conda-forge
bioconda
```

The `build` package is included in `environment.yml` so `python -m build` is available locally and in CI.

## CLI

Primary command:

```bash
sprite-mask run \
  --samples samples.tsv \
  --threshold 10 \
  --out results \
  --work work \
  --threads 4
```

Optional inputs:

```bash
sprite-mask run \
  --samples samples.tsv \
  --threshold 10 \
  --targets targets.bed \
  --mapq 20 \
  --reference ref.fa \
  --out results \
  --work work \
  --threads 4 \
  --keep-work
```

Recommended options:

- `--samples`: required TSV with sample metadata and alignment paths.
- `--threshold`: required integer depth threshold.
- `--out`: output directory; create if missing.
- `--work`: working directory; default to `<out>/work`.
- `--threads`: `mosdepth` decompression threads per sample.
- `--targets`: optional BED; if provided, clip each per-sample pass BED before counting and cohort intersection.
- `--mapq`: optional mapping-quality threshold passed to `mosdepth`.
- `--exclude-flag`: optional read flag exclusion value passed to `mosdepth`.
- `--reference`: required only when CRAM inputs need an explicit FASTA.
- `--strict-depth`: omit `mosdepth --fast-mode`.
- `--keep-work`: keep intermediates.
- `--force`: overwrite existing outputs.

Debug subcommands can be added later if helpful:

- `sprite-mask extract-pass`
- `sprite-mask collapse-samples`
- `sprite-mask collapse-populations`
- `sprite-mask summarize`

The first implementation should prioritize `sprite-mask run`.

## Input Contracts

### `samples.tsv`

Tab-delimited, with an optional header:

```text
sample_id  population  alignment
sample_1   popA        /path/sample_1.bam
sample_2   popA        /path/sample_2.bam
sample_3   popB        /path/sample_3.cram
```

Rules:

- `sample_id` values must be unique.
- `population` values are used as output column names in the population count BED.
- If the user does not care about populations, use `all` for every sample.
- `alignment` can be BAM or CRAM.

### `targets.bed`

Optional BED input. Only the first three columns are used:

```text
chrom  start  end
```

The tool should normalize this to a three-column working BED before using it.

## Test Dataset

Use a small regional subset of the 1000 Genomes high-coverage GRCh38 data for end-to-end tests.

Directory:

```text
tests/test_data/1000g_10sample_highcov_subset/
  samples.tsv
  sample_populations_and_sources.tsv
  samples.list
  targets.bed
  README.md
  bams/
    HG00096.bam
    HG00096.bam.bai
    ...
    NA18501.bam
    NA18501.bam.bai
  1000g_10samples_highcov.chr20_10000000-10100000.vcf.gz
  1000g_10samples_highcov.chr20_10000000-10100000.vcf.gz.tbi
  1000g_10samples_highcov.vcf.gz
  1000g_10samples_highcov.vcf.gz.tbi
```

Dataset description:

- Ten individuals from two populations: five GBR and five YRI.
- Source data are approximately 30x CRAMs from the 1000 Genomes/NYGC high-coverage collection.
- The source collection is GRCh38 and high coverage.
- The fixture BAMs are local regional BAM slices, optionally downsampled by the fixture-generation script to approximately 20x.
- The default region is `chr20:10000000-10100000` in GRCh38 coordinates.
- `targets.bed` contains the corresponding BED interval: `chr20 9999999 10100000`.

Primary files used by `sprite-mask` tests:

- `samples.tsv`: main metadata file for the CLI.
- `targets.bed`: three-column BED interval for the regional fixture.
- `bams/*.bam` and `bams/*.bam.bai`: one indexed local regional BAM per sample.

`samples.tsv` columns:

```text
sample_id  population  alignment
```

The `alignment` column points to the local regional BAM for each sample. Populations are `GBR` or `YRI`.

Additional fixture files:

- `sample_populations_and_sources.tsv`: provenance and source metadata for the fixture.
- `samples.list`: sample ID list.
- `1000g_10samples_highcov.chr20_10000000-10100000.vcf.gz`: indexed multi-sample VCF for the same region and the same ten samples.
- `1000g_10samples_highcov.vcf.gz`: symlink to the region-specific VCF.
- `1000g_10samples_highcov.vcf.gz.tbi`: symlink to the region-specific VCF index.

The VCF is included for fixture realism and possible future VCF-based tests. `sprite-mask` itself should primarily use `samples.tsv`, `targets.bed`, and the BAMs.

Verified fixture state:

- `samples.tsv` has 10 sample rows and 10 unique sample IDs.
- Population counts are exactly `GBR=5` and `YRI=5`.
- Every alignment path in `samples.tsv` exists.
- `bams/` contains 10 BAM files and 10 matching `.bam.bai` index files.
- BAM basenames match the sample IDs in `samples.tsv`.
- `samples.list` matches the sample order in `samples.tsv`.
- The VCF sample list matches the sample order in `samples.tsv`.
- The VCF contains records in the target interval and is readable with `bcftools`.
- The VCF header currently emits a `bcftools` warning that `AC` should be declared as `Number=A`; this does not affect the BAM-focused `sprite-mask` tests.
- `htsfile` recognizes the BAMs as BAM and the VCF as BGZF-compressed VCF.

Tooling note:

- Full end-to-end workflow validation requires `mosdepth` and `bedtools`, which should be available after creating the WSL `sprite` conda environment from `environment.yml`.
- The local macOS shell used for fixture verification had `bcftools` and `htsfile`, but not `mosdepth`, `bedtools`, or `samtools`.

## Output Contracts

Given `--threshold 10`, write:

```text
results/
  sample_sites_at_depth_10.tsv
  cohort.d10.sample_count_quantized.bed
  cohort.d10.population_count_quantized.bed
  cohort.d10.sites_by_sample_count.tsv
  cohort.d10.population_sites_by_count.tsv
  run_manifest.json
```

Working files:

```text
work/
  sample_1.d10.quantized.bed.gz
  sample_1.d10.mosdepth.summary.txt
  sample_1.d10.mosdepth.global.dist.txt
  sample_1.d10.pass.bed
  sample_1.d10.pass.targets.bed
  cohort.d10.multiinter.tsv
```

### `sample_sites_at_depth_10.tsv`

```text
sample_id  population  passing_sites
sample_1   popA        1234567
sample_2   popA        1200000
sample_3   popB        1180000
```

### `cohort.d10.sample_count_quantized.bed`

BED-like, sparse. Missing intervals mean zero passing samples.

```text
chrom  start  end  passing_samples
chr1   100    200  8
chr1   200    350  9
chr1   900    950  3
```

### `cohort.d10.population_count_quantized.bed`

BED-like, sparse. One count column per population.

```text
chrom  start  end  popA  popB  popC
chr1   100    200  4     3     1
chr1   200    350  5     3     1
chr1   900    950  2     1     0
```

### Summary TSVs

Sample-count summary:

```text
passing_samples  sites
1                1000
2                5000
3                9000
```

Population-count summary:

```text
population  passing_samples  sites
popA        0                1200
popA        1                5400
popA        2                9100
popB        0                3000
popB        1                7600
```

## Core Workflow

### 1. Parse and Validate Inputs

Module: `samples.py`

Functions:

```python
def read_samples(path: Path) -> list[Sample]:
    ...

def validate_samples(samples: list[Sample]) -> None:
    ...
```

Validation:

- Required columns are present or infer three columns when no header exists.
- Sample IDs are unique.
- Alignment paths exist.
- Threshold is a non-negative integer.
- Output paths are writable.
- `mosdepth` and `bedtools` are available on `PATH`.

Do not require `genome.txt`.

### 2. Normalize Targets

Module: `bedio.py`

Function:

```python
def normalize_targets_bed(targets_bed: Path, out_bed: Path) -> Path:
    ...
```

Behavior:

- Read `targets.bed` line by line.
- Skip blank and comment lines.
- Write only `chrom`, `start`, `end`.
- Validate `start < end`.
- Leave sorting to `bedtools sort`.

If no targets are supplied, skip this step.

### 3. Run Mosdepth Per Sample

Module: `mosdepth.py`

Function:

```python
def run_mosdepth(sample: Sample, config: RunConfig) -> MosdepthOutputs:
    ...
```

Command shape:

```bash
MOSDEPTH_Q0=FAIL \
MOSDEPTH_Q1=PASS \
mosdepth \
  --threads "${threads}" \
  --no-per-base \
  --quantize "0:${threshold}:" \
  "${prefix}" \
  "${alignment}"
```

Add `--fast-mode` unless `--strict-depth` is set.

Optional arguments:

- `--mapq <mapq>`
- `--flag <exclude_flag>`
- `--fasta <reference>` for CRAM/reference use

Expected outputs:

- `<prefix>.quantized.bed.gz`
- `<prefix>.mosdepth.summary.txt`
- `<prefix>.mosdepth.global.dist.txt`

Implementation notes:

- Use `subprocess.run(..., check=True)`.
- Pass `MOSDEPTH_Q0=FAIL` and `MOSDEPTH_Q1=PASS` through the environment.
- Capture stderr to a log file or manifest entry.
- Do not pipe `mosdepth` output through shell text tools.

### 4. Extract Per-Sample Pass BEDs

Module: `mosdepth.py` or `bedio.py`

Function:

```python
def extract_pass_intervals(quantized_bed_gz: Path, out_bed: Path) -> Path:
    ...
```

Behavior:

- Open the bgzipped BED with Python `gzip.open(..., "rt")`.
- Stream line by line.
- Keep rows where column 4 is `PASS`.
- Write only `chrom`, `start`, `end`.
- Do not load the file into memory.

Then sort and merge:

```python
def sort_and_merge_bed(in_bed: Path, out_bed: Path) -> Path:
    ...
```

Command shape:

```bash
bedtools sort -i input.bed | bedtools merge -i - > output.bed
```

Use subprocess piping in Python, not a shell string.

### 5. Clip Pass BEDs To Targets

Module: `bedtools.py`

Function:

```python
def intersect_sort_merge(a_bed: Path, b_bed: Path, out_bed: Path) -> Path:
    ...
```

Command shape:

```bash
bedtools intersect -a pass.bed -b targets.3col.bed |
bedtools sort -i - |
bedtools merge -i - > pass.targets.bed
```

Use target-clipped pass BEDs for both:

- Per-sample passing-site counts.
- `bedtools multiinter`.

If no targets are supplied, use the un-targeted pass BEDs.

### 6. Count Passing Sites Per Sample

Module: `bedio.py`

Function:

```python
def count_bed_sites(path: Path) -> int:
    total = 0
    for chrom, start, end in iter_bed3(path):
        total += end - start
    return total
```

Output:

```python
def write_sample_site_counts(samples: list[Sample], counts: dict[str, int], out_tsv: Path) -> None:
    ...
```

Write columns:

- `sample_id`
- `population`
- `passing_sites`

### 7. Run Bedtools Multiinter

Module: `bedtools.py`

Function:

```python
def run_multiinter(samples: list[Sample], pass_beds: list[Path], out_tsv: Path) -> Path:
    ...
```

Command shape:

```bash
bedtools multiinter \
  -header \
  -names sample_1 sample_2 sample_3 \
  -i sample_1.pass.bed sample_2.pass.bed sample_3.pass.bed \
  > cohort.d10.multiinter.tsv
```

Important:

- Do not pass `-empty`.
- Do not require or pass `-g`.
- Sparse output means absent intervals have zero passing samples.

Output columns:

```text
chrom  start  end  num  list  sample_1  sample_2  sample_3
```

### 8. Collapse By Total Sample Count

Module: `collapse.py`

Function:

```python
def collapse_sample_counts(multiinter_tsv: Path, out_bed: Path) -> None:
    ...
```

Algorithm:

1. Read header.
2. For each data row, parse `chrom`, `start`, `end`, and `num`.
3. Maintain a current interval `(chrom, start, end, count)`.
4. If the next row has the same `chrom`, starts at the current `end`, and has the same count, extend the current interval.
5. Otherwise flush the current interval and start a new one.
6. Write `chrom`, `start`, `end`, `passing_samples`.

This output is sparse because the input `multiinter` table is sparse.

### 9. Collapse By Population Count

Module: `collapse.py`

Function:

```python
def collapse_population_counts(samples: list[Sample], multiinter_tsv: Path, out_bed: Path) -> None:
    ...
```

Algorithm:

1. Build `sample_id -> population`.
2. Preserve population order from `samples.tsv`.
3. Read the `multiinter` header.
4. Validate that every sample column after column 5 maps to a population.
5. For each row, sum per-sample 0/1 indicators by population.
6. Maintain a current interval `(chrom, start, end, counts_tuple)`.
7. Merge adjacent rows only when the population-count tuple is identical.
8. Write a header:

```text
chrom  start  end  popA  popB  popC
```

This output is also sparse. Missing intervals mean zero samples pass in every population.

### 10. Write Summaries

Module: `summaries.py`

Functions:

```python
def summarize_sample_count_bed(sample_count_bed: Path, out_tsv: Path) -> None:
    ...

def summarize_population_count_bed(population_count_bed: Path, out_tsv: Path) -> None:
    ...
```

Sample-count summary:

- Read `cohort.d10.sample_count_quantized.bed`.
- Add `end - start` to bucket `passing_samples`.
- Write sorted count buckets.

Population-count summary:

- Read header from `cohort.d10.population_count_quantized.bed`.
- For each population column, add `end - start` to that population's bucket for the observed count.
- Write columns `population`, `passing_samples`, `sites`.

Because output is sparse, zero-count intervals are not included in these summaries.

### 11. Manifest And Logging

Module: `workflow.py`

Write `run_manifest.json`:

```json
{
  "threshold": 10,
  "targets": "targets.bed",
  "samples": 3,
  "threads": 4,
  "strict_depth": false,
  "outputs": {
    "sample_counts": "sample_sites_at_depth_10.tsv",
    "sample_count_bed": "cohort.d10.sample_count_quantized.bed",
    "population_count_bed": "cohort.d10.population_count_quantized.bed"
  },
  "commands": [
    {
      "tool": "mosdepth",
      "sample": "sample_1",
      "returncode": 0
    }
  ]
}
```

Keep the manifest simple at first. It should be enough to recover what was run.

## Implementation Order

### Phase 1: Package Skeleton

1. Create `environment.yml` for the WSL `sprite` conda environment.
2. Create `pyproject.toml`.
3. Create `src/sprite_mask/`.
4. Add `cli.py` with `argparse`.
5. Add `models.py` with dataclasses:

```python
@dataclass(frozen=True)
class Sample:
    sample_id: str
    population: str
    alignment: Path

@dataclass(frozen=True)
class RunConfig:
    threshold: int
    out_dir: Path
    work_dir: Path
    threads: int
    targets: Path | None
    mapq: int | None
    exclude_flag: int | None
    reference: Path | None
    strict_depth: bool
    keep_work: bool
    force: bool
```

### Phase 2: Pure Python File Logic

Implement and test:

1. `read_samples`.
2. `normalize_targets_bed`.
3. `extract_pass_intervals`.
4. `count_bed_sites`.
5. `collapse_sample_counts`.
6. `collapse_population_counts`.
7. Summary writers.

This phase should use tiny fixture files and no external tools.

### Phase 3: Subprocess Wrappers

Implement:

1. `check_tool("mosdepth")`.
2. `check_tool("bedtools")`.
3. `run_command`.
4. `run_mosdepth`.
5. `sort_and_merge_bed`.
6. `intersect_sort_merge`.
7. `run_multiinter`.

Tests can mock `subprocess.run` for command construction.

### Phase 4: End-To-End Workflow

Implement `run_workflow(config)`:

1. Parse samples.
2. Normalize targets if provided.
3. For each sample:
   - run `mosdepth`
   - extract pass intervals
   - sort/merge pass intervals
   - clip to targets if provided
   - count passing sites
4. Write per-sample site counts.
5. Run `bedtools multiinter`.
6. Collapse sample counts.
7. Collapse population counts.
8. Write summaries.
9. Write manifest.
10. Clean work directory unless `--keep-work`.

### Phase 5: Minimal Integration Test

Use two levels of integration testing.

First, keep a tiny synthetic BED fixture for fast interval logic tests:

```text
sample_1.pass.bed
chr1  0   10
chr1  20  30

sample_2.pass.bed
chr1  5   25

sample_3.pass.bed
chr1  20  30
```

Expected sparse `multiinter` shape:

```text
chr1  0   5   1
chr1  5   10  2
chr1  10  20  1
chr1  20  25  3
chr1  25  30  2
```

Expected collapsed sample-count BED:

```text
chr1  0   5   1
chr1  5   10  2
chr1  10  20  1
chr1  20  25  3
chr1  25  30  2
```

With populations:

```text
sample_1 popA
sample_2 popA
sample_3 popB
```

Expected population-count BED:

```text
chrom  start  end  popA  popB
chr1   0      5    1     0
chr1   5      10   2     0
chr1   10     20   1     0
chr1   20     25   2     1
chr1   25     30   1     1
```

Second, add an end-to-end fixture test using:

```text
tests/test_data/1000g_10sample_highcov_subset/
```

Test command:

```bash
sprite-mask run \
  --samples tests/test_data/1000g_10sample_highcov_subset/samples.tsv \
  --targets tests/test_data/1000g_10sample_highcov_subset/targets.bed \
  --threshold 10 \
  --out tmp/1000g_fixture_results \
  --threads 2 \
  --keep-work
```

Expected checks:

- The command exits successfully.
- `sample_sites_at_depth_10.tsv` contains ten samples.
- The populations represented are exactly `GBR` and `YRI`.
- `cohort.d10.sample_count_quantized.bed` exists and is non-empty.
- `cohort.d10.population_count_quantized.bed` exists and has columns `chrom`, `start`, `end`, `GBR`, and `YRI`.
- All output intervals fall within `targets.bed`.
- Passing-sample counts are between 1 and 10.
- Population count columns are between 0 and 5.
- The population-count BED and sample-count BED are consistent: after splitting both outputs at their combined breakpoints, `GBR + YRI` equals `passing_samples` for every overlapping sub-interval.
- No CRAM reference lookup is needed because the fixture uses local indexed BAMs.

### Phase 6: Cross-Platform Build CI

Add CI that runs on Linux and macOS for every pull request and release tag.

Required jobs:

1. `test-linux`: create `sprite`, run linting, type checks, unit tests, and integration tests.
2. `test-macos`: same as Linux.
3. `build-linux`: build sdist and wheel after tests pass.
4. `build-macos`: build sdist and wheel after tests pass.

If the package remains pure Python, the wheel artifact from either OS should be equivalent. Still run the build step on both OSes because it catches packaging and environment problems early.

## Testing Plan

Unit tests:

- `read_samples` accepts header and no-header files.
- `read_samples` rejects duplicate sample IDs.
- `extract_pass_intervals` keeps only `PASS`.
- `count_bed_sites` sums `end - start`.
- `collapse_sample_counts` merges adjacent rows with the same count.
- `collapse_sample_counts` does not merge non-adjacent rows with the same count.
- `collapse_population_counts` preserves population order.
- `collapse_population_counts` merges only identical population-count vectors.
- Summary functions write sorted buckets.

Command-wrapper tests:

- `run_mosdepth` builds the expected command and environment.
- `run_mosdepth` includes `--fast-mode` by default.
- `run_mosdepth` omits `--fast-mode` with `--strict-depth`.
- `run_multiinter` does not include `-empty` or `-g`.
- `run_multiinter` passes sample names in `samples.tsv` order.

Integration tests:

- If `bedtools` is available, run `multiinter` on tiny pass BED fixtures.
- Compare collapsed outputs to expected files.
- Mark these tests skipped when `bedtools` is missing.
- If `mosdepth` and `bedtools` are available, run the full workflow on `tests/test_data/1000g_10sample_highcov_subset`.
- Keep the 1000G fixture test marked as an integration test so it can be skipped in fast unit-only runs.
- Run integration tests on both Linux and macOS CI jobs.

Build tests:

- Run `python -m build` on Linux and macOS.
- Install the built wheel into a fresh `sprite` environment when practical.
- Run `sprite-mask --help` and `sprite-mask run --help` from the installed wheel.

## Error Handling

Keep errors direct and actionable:

- Missing `mosdepth`: `mosdepth not found on PATH`.
- Missing `bedtools`: `bedtools not found on PATH`.
- Empty `samples.tsv`: `samples.tsv contains no samples`.
- Duplicate sample: `duplicate sample_id: sample_1`.
- Missing alignment: `alignment does not exist for sample_1: /path/file.bam`.
- Bad BED interval: `invalid BED interval at targets.bed:12`.
- Missing `PASS` output: `mosdepth did not produce expected quantized BED for sample_1`.

## Performance Choices

- Stream all TSV/BED parsing.
- Keep `mosdepth --no-per-base`.
- Use `mosdepth --fast-mode` by default.
- Avoid shell text processing.
- Avoid materializing zero-pass intervals.
- Keep raw `mosdepth` quantized files for reproducibility unless cleanup is explicitly requested later.

## First Milestone

The first usable version is complete when this command works:

```bash
sprite-mask run \
  --samples samples.tsv \
  --threshold 10 \
  --out results \
  --threads 4
```

And writes:

- `sample_sites_at_depth_10.tsv`
- `cohort.d10.sample_count_quantized.bed`
- `cohort.d10.population_count_quantized.bed`
- `cohort.d10.sites_by_sample_count.tsv`
- `cohort.d10.population_sites_by_count.tsv`
- `run_manifest.json`
