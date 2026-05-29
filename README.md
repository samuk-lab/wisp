``sprite``<img src="https://raw.githubusercontent.com/samuk-lab/sprite/master/docs/images/sprite_logo.png" align="right" width="20%">
====================

`sprite` is a CLI for building population count masks from BAM/CRAM/VCF files. It is a companion
utility for [pixy](https://github.com/ksamuk/pixy/), but can also be used on its own!

population count masks can be used along with a variants-only VCF to properly compute the denominators
of pi, dxy, Watterson's theta, and Tajima's D.

The package is pre-release alpha pending validation, but will be distributed on bioconda as `sprite-mask` (just `sprite` was taken!).

## Install (pre-release instructions)

### 1. Make a clean environment
mamba create -n sprite python=3.11 pip git -c conda-forge
mamba activate sprite

### 2. Install command-line dependencies
mamba install -c conda-forge samtools bcftools htslib mosdepth

### 3. Install sprite directly from GitHub
python -m pip install "git+https://github.com/samuk-lab/sprite.git"

## Install (when the bioconda version becomes available):

Create the conda environment from the repository root:

```bash
mamba env create -f environment.yml
conda activate sprite
python -m pip install -e ".[dev]"
```

If the environment already exists:

```bash
mamba env update -n sprite -f environment.yml
conda activate sprite
python -m pip install -e ".[dev]"
```

## Run

```bash
sprite from-alignments \
  --samples tests/test_data/1000g_5sample_chr20_smoke/samples.tsv \
  --threshold 10 \
  --targets tests/test_data/1000g_5sample_chr20_smoke/targets.bed \
  --out results \
  --work work \
  --threads 4 \
  --jobs 2 \
  --keep-work
```

The sample table is tab-delimited and may have a header:

```text
sample_id  population  alignment
sample_1   popA        /path/sample_1.bam
sample_2   popB        /path/sample_2.cram
```

If BAM/CRAM read group sample names are present, they must match the
corresponding `sample_id`.

To build the same output from an all-sites VCF:

```bash
sprite from-vcf \
  --all-sites-vcf validation/1000g_20sample_highcov_4chrom_subset/1000g_20samples_highcov_4chroms.all_sites.bam_call.trim_alt.vcf.gz \
  --popfile validation/1000g_20sample_highcov_4chrom_subset/sample_populations.tsv \
  --threshold 10 \
  --out results \
  --work work
```

The popfile is tab-delimited and may have a header:

```text
sample_id  population
sample_1   popA
sample_2   popB
```

VCF mode requires every popfile sample to be present in the VCF. VCF sample
columns that are absent from the popfile produce a warning and are ignored.
Records may carry VCF `FILTER` values; the input is assumed to have already been
filtered as desired.
When duplicate records share a `CHROM:POS`, a sample passes that base if any duplicate
record has `FORMAT/DP >= --threshold`. Duplicate records must be contiguous, as in
a coordinate-sorted VCF.

Regardless of mode, sprite outputs two main files:

- `sprite.bed.gz`
- `sprite.bed.gz.tbi`

Set `--output-prefix` to choose a different output filename prefix; `.bed.gz`
is appended to the prefix.

The final population BED output is sparse: intervals absent from the file are interpreted as
zero passing samples per population. The `*.bed.gz` file
includes comment-prefixed headers. The first header line contains JSON metadata
with the output format, threshold, columns, and explicit population-column
mappings/sample counts where applicable; the second header line is the column
header, for example `#chrom	start	end	GBR	YRI`. The `*.bed.gz`
file is sorted, bgzip-compressed, and tabix-indexed.

## Test Data

Fixture-generation scripts using the Hg1000 dataset live in `tests/test_data/scripts/`:

```bash
bash tests/test_data/scripts/download_1000g_10sample_chr20_highcov_fixture.sh
bash tests/test_data/scripts/download_1000g_20sample_4chrom_highcov_fixture.sh
```

## Documentation

The project documentation lives in `docs/`.
