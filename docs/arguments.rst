Arguments
*********

All arguments are listed below. ``sprite --help`` shows the same information.

Commands
========

``sprite`` has two subcommands:

**from-alignments**
    Build a population count mask from BAM/CRAM files via ``mosdepth``.

**from-vcf**
    Build a population count mask from a prefiltered all-sites VCF with
    per-sample ``FORMAT/DP`` values.

Required for every run
======================

**--min-dp INTEGER**
    Minimum depth for a sample to pass a site. Must be non-negative.
    ``--min-dp 0`` requires ``--mask`` because no genome-wide coordinate
    file is supplied.

**--out PATH**
    Output directory for the final ``sprite.bed.gz`` and tabix index.

Input-specific arguments
========================

**--samples PATH**
    Sample metadata TSV for BAM/CRAM mode. Must provide sample ID,
    population, and alignment path columns. See :doc:`inputs`.

**--all-sites-vcf PATH**
    All-sites VCF for VCF mode. Must include per-sample ``FORMAT/DP``
    values for every sample in ``--popfile``.

**--popfile PATH**
    Sample/population TSV for VCF mode. Required with ``--all-sites-vcf``.

Shared optional arguments
=========================

**--work PATH**
    Working directory for intermediate files. Defaults to ``<out>/work``.

**--mask PATH**
    Restrict output to sites within this BED. In BAM/CRAM mode, the mask
    is normalized to three columns, sorted, and merged before sample pass
    intervals are clipped. In VCF mode, only bases within the mask
    intervals are emitted.

**--output-prefix TEXT**
    Output filename stem within ``--out``. Defaults to ``sprite``,
    producing ``sprite.bed.gz`` and ``sprite.bed.gz.tbi``. ``.bed.gz``
    is always appended.

**--keep-work**
    Keep intermediate files in the working directory. By default, work
    files are removed after a successful run.

**--force**
    Overwrite existing final outputs. Without this flag, ``sprite`` refuses
    to replace ``sprite.bed.gz`` or ``sprite.bed.gz.tbi``.

**--version**
    Print the installed ``sprite`` version and exit.

**--help**
    Print the full help message and exit.

BAM/CRAM mode options
=====================

**--threads INTEGER**
    ``mosdepth`` threads per sample. Defaults to ``1``.

**--jobs INTEGER**
    Samples to process concurrently. Defaults to ``1``. Total parallelism
    is ``--jobs × --threads``.

**--min-mapq INTEGER**
    Minimum read mapping quality.

**--max-dp INTEGER**
    Maximum depth to pass a site.

**--exclude-flag INTEGER**
    SAM FLAG bits to exclude reads.

**--reference PATH**
    FASTA reference for CRAM inputs.

**--strict-depth**
    Use precise per-base depth counting. Slower, but avoids the
    approximations made by ``mosdepth --fast-mode``.

Examples
========

BAM/CRAM mode:

.. code-block:: console

   sprite from-alignments \
     --samples tests/test_data/1000g_5sample_chr20_smoke/samples.tsv \
     --min-dp 10 \
     --mask tests/test_data/1000g_5sample_chr20_smoke/targets.bed \
     --out results \
     --work work \
     --threads 4 \
     --jobs 2 \
     --keep-work

All-sites VCF mode:

.. code-block:: console

   sprite from-vcf \
     --all-sites-vcf validation/cohort.all_sites.vcf.gz \
     --popfile validation/sample_populations.tsv \
     --min-dp 10 \
     --mask validation/targets.bed \
     --out results \
     --work work
