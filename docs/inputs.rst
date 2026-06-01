Input files
***********

Sample table for BAM/CRAM mode
===============================

``--samples`` points to a tab-delimited table with three required fields:
sample ID, population, and alignment path.

The table accepts an optional header:

.. code-block:: text

   sample_id  population  alignment
   HG00096    GBR         bams/HG00096.bam
   NA18486    YRI         bams/NA18486.bam

Accepted header aliases:

* Sample ID: ``sample_id``, ``sample``, or ``id``
* Population: ``population`` or ``pop``
* Alignment: ``alignment``, ``bam_or_cram``, ``bam``, or ``cram``

Without a recognized header, columns are read in order:
``sample_id``, ``population``, ``alignment``.

Alignment paths may be absolute or relative. For relative paths, ``sprite``
first checks the current value, then the path relative to the sample table's
parent directory.

When read group sample names are present in a BAM/CRAM header, each ``@RG``
``SM`` value must match the ``sample_id`` for that row.

Popfile for all-sites VCF mode
===============================

``--popfile`` points to a tab-delimited table with sample ID and population
columns:

.. code-block:: text

   sample_id  population
   HG00096    GBR
   NA18486    YRI

The same header aliases apply. Without a recognized header, columns are read
as ``sample_id`` then ``population``.

Sample and population labels
============================

Sample IDs and population names must be non-empty and must not contain
whitespace or path separators. Sample IDs must be unique.

All-sites VCF
=============

``--all-sites-vcf`` must point to a VCF with a ``#CHROM`` header and sample
columns. Every sample ID in ``--popfile`` must appear in the VCF. VCF samples
absent from ``--popfile`` produce a warning and are ignored.

For each record, ``sprite`` reads the ``DP`` field from the sample's
``FORMAT`` value. Missing DP values do not pass. Non-integer DP values are
rejected. Records may carry any ``FILTER`` value; the file is assumed to have
been filtered as desired. Duplicate ``CHROM:POS`` records are merged with OR
semantics per sample, but duplicates must be contiguous, as in a
coordinate-sorted VCF.

Variants-only VCF for BAM/CRAM mode
===================================

``--variants-vcf`` is optional in BAM/CRAM mode. It should point to a
coordinate-sorted variants-only VCF for the same samples and reference
coordinate system as the alignments. When sample columns are present, every
sample ID in ``--samples`` must appear in the VCF; extra VCF samples are
ignored for threshold estimation.

``sprite`` uses this VCF in two ways.

Threshold estimation
--------------------

If ``--min-dp`` or ``--max-dp`` is omitted, ``sprite`` estimates it from
positive per-sample ``FORMAT/DP`` values at variant records. ``--min-dp``
uses the smallest observed positive DP and ``--max-dp`` uses the largest
observed positive DP. If ``--min-mapq`` is omitted, ``sprite`` estimates it
from the smallest ``INFO/MQ`` value, rounded down to an integer.

Any threshold supplied manually on the command line takes precedence over
the VCF-derived estimate. If ``--min-dp`` is omitted and the VCF has no usable
``FORMAT/DP`` values, the run is rejected.

Variant exclusions
------------------

The same VCF is scanned for variant classes that should not contribute
callable single-base denominators: indels, symbolic structural variants,
breakends, and multi-nucleotide polymorphisms. SNP-only records are retained.

Exclusion spans are emitted in BED coordinates. For symbolic structural
variants, ``INFO/END`` is preferred when present; otherwise ``SVLEN`` is used
when available. For ordinary sequence alleles, the reference allele length is
used. The resulting intervals are sorted, merged, and subtracted from every
sample pass BED before population counts are built.

Mask BED
========

``--mask`` accepts a BED file with at least three columns. Blank lines,
comments, and a ``chrom start end`` header are ignored. Coordinates must be
non-negative, 0-based, half-open intervals where ``start < end``.

In BAM/CRAM mode, the mask is sorted and merged with ``bedtools`` before use.
In VCF mode, overlapping intervals are merged internally before bases are
tested for containment.
