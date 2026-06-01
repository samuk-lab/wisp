Examples
********

BAM/CRAM smoke run
==================

This run uses the small fixture bundled with the test suite:

.. code-block:: console

   sprite from-alignments \
     --samples tests/test_data/1000g_5sample_chr20_smoke/samples.tsv \
     --min-dp 10 \
     --mask tests/test_data/1000g_5sample_chr20_smoke/targets.bed \
     --out results/chr20_smoke \
     --work work/chr20_smoke \
     --threads 2 \
     --jobs 2

Expected outputs:

.. code-block:: text

   results/chr20_smoke/sprite.bed.gz
   results/chr20_smoke/sprite.bed.gz.tbi

Keep intermediate files
=======================

Add ``--keep-work`` to inspect the mosdepth outputs, sample pass BEDs, and
``bedtools multiinter`` results:

.. code-block:: console

   sprite from-alignments \
     --samples tests/test_data/1000g_5sample_chr20_smoke/samples.tsv \
     --min-dp 10 \
     --mask tests/test_data/1000g_5sample_chr20_smoke/targets.bed \
     --out results/chr20_debug \
     --work work/chr20_debug \
     --keep-work

BAM/CRAM with a variants-only VCF
=================================

Use ``--variants-vcf`` when you have a variants-only VCF from the same callset
and want ``sprite`` to estimate omitted alignment thresholds and mask
non-SNP variant spans:

.. code-block:: console

   sprite from-alignments \
     --samples tests/test_data/1000g_5sample_chr20_smoke/samples.tsv \
     --variants-vcf validation/cohort.variants.vcf.gz \
     --mask tests/test_data/1000g_5sample_chr20_smoke/targets.bed \
     --out results/chr20_variants_vcf \
     --work work/chr20_variants_vcf \
     --threads 2 \
     --jobs 2

Manual threshold flags override VCF-derived estimates:

.. code-block:: console

   sprite from-alignments \
     --samples samples.tsv \
     --variants-vcf cohort.variants.vcf.gz \
     --min-dp 8 \
     --max-dp 80 \
     --min-mapq 30 \
     --out results/manual_thresholds

All-sites VCF run
=================

When you have a prefiltered all-sites VCF with per-sample DP values:

.. code-block:: console

   sprite from-vcf \
     --all-sites-vcf validation/1000g_20sample_highcov_4chrom_subset/1000g_20samples_highcov_4chroms.all_sites.bam_call.trim_alt.vcf.gz \
     --popfile validation/1000g_20sample_highcov_4chrom_subset/sample_populations.tsv \
     --min-dp 10 \
     --mask validation/1000g_20sample_highcov_4chrom_subset/targets.bed \
     --out results/vcf_mode \
     --work work/vcf_mode

Query the output
================

The final BED is tabix-indexed, so you can slice it directly:

.. code-block:: console

   tabix results/chr20_smoke/sprite.bed.gz chr20:10000000-10010000

Note that BED coordinates are 0-based and half-open, while tabix region
strings are 1-based inclusive.

Overwrite existing outputs
==========================

``sprite`` refuses to overwrite final outputs unless ``--force`` is passed:

.. code-block:: console

   sprite from-alignments \
     --samples tests/test_data/1000g_5sample_chr20_smoke/samples.tsv \
     --min-dp 10 \
     --mask tests/test_data/1000g_5sample_chr20_smoke/targets.bed \
     --out results/chr20_smoke \
     --force
