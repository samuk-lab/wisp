About
*****

``sprite`` creates sparse depth masks for population genomic workflows.
Instead of reporting a single cohort-wide pass/fail mask, it reports the number
of passing samples per population across genomic intervals.

Why population-count masks?
===========================

Downstream analyses often need to know whether enough individuals in each
population have usable depth at a site or region. A population-count BED keeps
that information explicit:

.. code-block:: text

   #chrom  start  end  GBR  YRI
   chr20   9999999  10001042  7  4
   chr20   10001042 10001881 10  6

Each population column contains a count of samples passing the configured depth
threshold over the interval. Adjacent bases with the same population counts are
collapsed into one BED interval.

Input modes
===========

BAM/CRAM mode
-------------

In alignment mode, ``sprite`` runs ``mosdepth`` once per sample with a two-bin
quantization: below threshold and passing threshold. It extracts the passing
intervals, optionally clips them to a target BED, intersects all sample pass BEDs
with ``bedtools multiinter``, and collapses per-sample indicators into
population counts.

All-sites VCF mode
------------------

In VCF mode, ``sprite`` reads ``FORMAT/DP`` values directly from an all-sites
VCF. Each requested sample passes a base when its DP value is greater than or
equal to ``--threshold``. Duplicate records at the same ``CHROM:POS`` are merged
with OR semantics for each sample: if any duplicate record passes, the sample
passes that base. Duplicate records must be contiguous, as in a
coordinate-sorted VCF.

Sparse output
=============

``sprite`` omits intervals where all population counts are zero. This is why
the output is described as sparse. Consumers should treat missing intervals as
zero passing samples per population, not as unknown or missing output.
