About
*****

``sprite`` creates population count masks for population genomic workflows.
Where a conventional depth mask gives a single cohort-wide pass/fail per
site, ``sprite`` reports how many samples in each population clear the depth
threshold.

Why population count masks?
===========================

Downstream analyses often need to know whether enough individuals in each
population have usable depth at a site or region. A population count mask
makes that explicit:

.. code-block:: text

   #chrom  start  end  GBR  YRI
   chr20   9999999  10001042  7  4
   chr20   10001042 10001881 10  6

Each population column is a count of samples passing the configured depth
threshold over that interval. Adjacent bases with the same counts are
collapsed into a single BED interval.

Input modes
===========

BAM/CRAM mode
-------------

In alignment mode, ``sprite`` runs ``mosdepth`` once per sample using a
two-bin quantization: below threshold and at-or-above threshold. It extracts
the passing intervals, optionally clips them to a mask BED, intersects all
sample pass BEDs with ``bedtools multiinter``, and assembles them into a
population count mask.

All-sites VCF mode
------------------

In VCF mode, ``sprite`` reads ``FORMAT/DP`` values directly from an all-sites
VCF. A sample passes a base when its DP value is greater than or equal to
``--min-dp``. Duplicate records at the same ``CHROM:POS`` are merged with OR
semantics per sample: if any duplicate passes, the sample passes that base.
Duplicate records must be contiguous, as in a coordinate-sorted VCF.

Sparse output
=============

``sprite`` omits intervals where all population counts are zero. Consumers
should treat missing intervals as zero passing samples per population, not as
unknown or skipped.
