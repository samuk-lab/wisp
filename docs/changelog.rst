Changelog
*********

0.1.0
=====

Initial documented release.

Highlights:

* Build sparse population-count BEDs from BAM/CRAM cohorts.
* Use an optional variants-only VCF in BAM/CRAM mode to estimate omitted
  depth/MAPQ thresholds and exclude indel, structural-variant, breakend, and
  multi-nucleotide polymorphism spans.
* Build the same output from prefiltered all-sites VCF ``FORMAT/DP`` values.
* Write bgzipped and tabix-indexed ``sprite.bed.gz`` output.
* Include JSON metadata and population column headers in the output BED.
