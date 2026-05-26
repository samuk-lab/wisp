# 1000 Genomes 10-sample high-coverage test subset

This directory contains a small, regional 1000 Genomes high-coverage test dataset.

## Samples

- 5 GBR samples: HG00096, HG00097, HG00099, HG00100, HG00101
- 5 YRI samples: NA18486, NA18489, NA18498, NA18499, NA18501

## Region

`chr20:10000000-10100000`

Coordinates are GRCh38 / chr-prefixed 1000 Genomes high-coverage coordinates.

## Coverage

The source CRAMs are NYGC high-coverage, approximately 30x whole-genome data.
This script uses DOWNSAMPLE_FRAC=0.67, so the regional BAMs are
approximately 30x * 0.67. Set DOWNSAMPLE_FRAC=1 to keep full depth.

## Files

- `samples.tsv`: sprite-mask sample metadata
- `sample_populations_and_sources.tsv`: sample/population/source metadata
- `samples.list`: sample list used for VCF subsetting
- `targets.bed`: BED interval for the selected region
- `bams/*.bam`: regional BAMs for each sample
- `1000g_10samples_highcov.vcf.gz`: multi-sample VCF for all 10 samples in the selected region

## Note

The source alignments are CRAMs, but the local test files are BAMs so downstream
tests do not need to handle CRAM reference lookup.
