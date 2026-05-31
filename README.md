`sprite`<img src="https://raw.githubusercontent.com/samuk-lab/sprite/master/docs/images/sprite_logo.png" align="right" width="20%">
====================

`sprite` builds population count masks from BAM/CRAM alignments or all-sites VCFs. It is a companion tool for [pixy](https://github.com/ksamuk/pixy/), though it works equally well on its own.

Population count masks let you correctly compute the denominators of π, d<sub>xy</sub>, Watterson's θ, and Tajima's D when working from a variants-only VCF — callable sites are counted per population rather than collapsed into a single cohort-wide pass/fail.

> **Note:** `sprite` is pre-release (alpha) and pending validation. It will be distributed on Bioconda as `sprite-mask` once stable.

## Installation

```bash
mamba create -n sprite python=3.11 pip git -c conda-forge
mamba activate sprite
mamba install -c conda-forge samtools bcftools htslib mosdepth
python -m pip install "git+https://github.com/samuk-lab/sprite.git"
```

## Usage

### From BAM/CRAM alignments

`sprite` runs [mosdepth](https://github.com/brentp/mosdepth) on each sample and collapses per-sample pass intervals into a population count mask:

```bash
sprite from-alignments \
  --samples samples.tsv \
  --min-dp 10 \
  --out results \
  --threads 4 \
  --jobs 2
```

The sample table is tab-delimited with an optional header:

```text
sample_id  population  alignment
sample_1   popA        /path/sample_1.bam
sample_2   popB        /path/sample_2.cram
```

If BAM/CRAM read groups include sample names, they must match the corresponding `sample_id`.

### From an all-sites VCF

```bash
sprite from-vcf \
  --all-sites-vcf all_sites.vcf.gz \
  --popfile populations.tsv \
  --min-dp 10 \
  --out results
```

The population file is tab-delimited with an optional header:

```text
sample_id  population
sample_1   popA
sample_2   popB
```

A few things to know about VCF mode:

- Every sample in the population file must appear in the VCF. VCF samples absent from the population file are ignored (with a warning).
- Records may carry any `FILTER` value — the input is assumed to have been filtered as desired before running `sprite`.
- At duplicate `CHROM:POS` records, a sample passes a site if any duplicate has `FORMAT/DP ≥ --min-dp`. Duplicates must be contiguous, as in a coordinate-sorted VCF.
- By default, all record types are used (SNPs, indels, symbolic alleles, invariant sites). Pass `--snps-only` to exclude indel sites while retaining invariant sites.

## Output

`sprite` writes two files to `--out`:

| File | Description |
|---|---|
| `sprite.bed.gz` | bgzip-compressed, tabix-indexed population count mask |
| `sprite.bed.gz.tbi` | tabix index |

Use `--output-prefix` to change the filename stem; `.bed.gz` is always appended.

The population count mask is sparse: intervals absent from the output carry zero passing samples across all populations. Two comment-prefixed header lines precede the data:

1. A JSON metadata line recording run parameters, output format, column names, and per-population sample counts.
2. A column header line, e.g. `#chrom	start	end	GBR	YRI`.

## Test data

Fixture-generation scripts using 1000 Genomes data live in `tests/test_data/scripts/`:

```bash
bash tests/test_data/scripts/download_1000g_10sample_chr20_highcov_fixture.sh
bash tests/test_data/scripts/download_1000g_20sample_4chrom_highcov_fixture.sh
```

## Documentation

Full documentation lives in `docs/`.
