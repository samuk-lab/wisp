Arguments
*********

Below is the full list of arguments accepted by ``sprite``. You can also see
the canonical help text with ``sprite --help``.

Commands
========

``sprite`` has two run commands:

**from-alignments**
    Build a mask from BAM/CRAM files via ``mosdepth``.

**from-vcf**
    Build a mask from a prefiltered all-sites VCF with per-sample
    ``FORMAT/DP`` values.

Required for every run
======================

**--threshold INTEGER**
    Minimum depth required for a sample to pass an interval or base. The value
    must be non-negative. In BAM/CRAM mode, ``--threshold 0`` requires
    ``--targets`` because no genome-wide coordinate file is supplied.

**--out PATH**
    Output directory for the final ``cohort.sprite.bed.gz`` and tabix index.

Required input mode
===================

Use the command for the input mode you want.

**--samples PATH**
    Sample metadata TSV for BAM/CRAM mode. The table must provide sample ID,
    population, and alignment path columns. See :doc:`inputs`.

**--all-sites-vcf PATH**
    Prefiltered all-sites VCF for VCF mode. The VCF must include per-sample
    ``FORMAT/DP`` values for every sample in ``--popfile``.

**--popfile PATH**
    Sample/population TSV for VCF mode. Required with ``--all-sites-vcf`` and
    not valid with ``--samples``.

Shared optional arguments
=========================

**--work PATH**
    Working directory for intermediate files. Defaults to ``<out>/work``.

**--targets PATH**
    Optional BED regions to include. In BAM/CRAM mode, targets are normalized to
    three columns, sorted, and merged before sample pass intervals are clipped.
    In VCF mode, only bases contained by the target intervals are emitted.

**--keep-work**
    Keep intermediate files in the working directory. By default, generated
    work files are removed after a successful run.

**--force**
    Overwrite existing final outputs. Without this flag, ``sprite`` refuses to
    replace ``cohort.sprite.bed.gz`` or ``cohort.sprite.bed.gz.tbi``.

**--version**
    Print the installed ``sprite`` version and exit.

**--help**
    Print the full help message and exit.

BAM/CRAM mode options
=====================

**--threads INTEGER**
    Number of ``mosdepth`` threads per sample. Defaults to ``1``.

**--jobs INTEGER**
    Number of samples to process concurrently. Defaults to ``1``.

**--mapq INTEGER**
    Pass a ``mosdepth --mapq`` mapping-quality filter.

**--exclude-flag INTEGER**
    Pass a ``mosdepth --flag`` read flag exclusion value.

**--reference PATH**
    FASTA reference passed to ``mosdepth --fasta`` for CRAM inputs.

**--strict-depth**
    Omit ``mosdepth --fast-mode``. By default, ``sprite`` enables fast mode.

Examples
========

BAM/CRAM mode:

.. code-block:: console

   sprite from-alignments \
     --samples tests/test_data/1000g_5sample_chr20_smoke/samples.tsv \
     --threshold 10 \
     --targets tests/test_data/1000g_5sample_chr20_smoke/targets.bed \
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
     --threshold 10 \
     --targets validation/targets.bed \
     --out results \
     --work work
