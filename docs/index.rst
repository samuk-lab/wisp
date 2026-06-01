.. sprite documentation master file.

.. raw:: html

   <div align="center"><h1>sprite 0.1.0</h1></div>

.. image:: images/sprite_logo.png
   :width: 200
   :align: center

What is sprite?
===============

``sprite`` is a command line tool for building population count masks: sparse,
population-level summaries of how many individuals have sufficient depth and
mapping quality to call genotypes at sites across the genome.
These masks can be produced from BAM/CRAM alignments or an all-sites VCF.

.. list-table::
   :header-rows: 1

   * - chrom
     - start
     - end
     - GBR
     - YRI
   * - chr20
     - 9999999
     - 10001042
     - 7
     - 4
   * - chr20
     - 10001042
     - 10001881
     - 10
     - 6
   * - chr20
     - 10500000
     - 10501200
     - 3
     - 0

Each row records how many samples in each population cleared the depth
threshold over that interval. Intervals where all counts are zero are omitted.

Population count masks let you correctly compute the denominators of π, d\
:sub:`xy`, Watterson's θ, and Tajima's D when working from a variants-only
VCF — callable sites are counted per population rather than collapsed into a
single cohort-wide pass/fail.

``sprite`` is designed for use with `pixy <https://pixy.readthedocs.io>`_,
but works equally well on its own. Source code is on
`GitHub <https://github.com/samuk-lab/sprite>`_.

The tool produces the same ``sprite.bed.gz`` output from two input modes:

* BAM/CRAM alignments, using ``mosdepth`` to quantize each sample and
  ``bedtools multiinter`` to combine samples. In this mode, an optional
  variants-only VCF can estimate omitted thresholds and exclude non-SNP
  variant spans from the final sparse mask.
* A prefiltered all-sites VCF, using per-sample ``FORMAT/DP`` values directly.

The output is bgzip-compressed and tabix-indexed. Large cohorts and large
target regions stay compact because absent intervals are interpreted as zero
passing samples in every population.

.. toctree::
   :caption: Documentation
   :maxdepth: 2

   about
   installation
   arguments
   inputs
   examples
   output

.. toctree::
   :caption: Reference
   :maxdepth: 2

   development
   api
   changelog
