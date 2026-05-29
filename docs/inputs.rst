Input Files
***********

Sample table for BAM/CRAM mode
==============================

``--samples`` points to a tab-delimited table with three required fields:
sample ID, population, and alignment path.

The table may include a header:

.. code-block:: text

   sample_id  population  alignment
   HG00096    GBR         bams/HG00096.bam
   NA18486    YRI         bams/NA18486.bam

Header aliases are accepted:

* Sample ID: ``sample_id``, ``sample``, or ``id``
* Population: ``population`` or ``pop``
* Alignment: ``alignment``, ``bam_or_cram``, ``bam``, or ``cram``

Without a recognized header, columns are interpreted in this order:
``sample_id``, ``population``, ``alignment``.

Alignment paths may be absolute or relative. For relative paths, ``sprite``
first checks the current value, then the path relative to the sample table's
parent directory.

Popfile for all-sites VCF mode
==============================

``--popfile`` points to a tab-delimited table with sample ID and population
columns:

.. code-block:: text

   sample_id  population
   HG00096    GBR
   NA18486    YRI

The same sample and population header aliases are accepted. Without a
recognized header, columns are interpreted as ``sample_id`` then ``population``.

Sample and population labels
============================

Sample IDs and population names must be non-empty. They must not contain
whitespace or path separators. Sample IDs must be unique.

All-sites VCF
=============

``--all-sites-vcf`` must point to a VCF with a ``#CHROM`` header and sample
columns. Every sample in ``--popfile`` must be present in the VCF.

For each record, ``sprite`` reads the ``DP`` field from the sample's
``FORMAT`` value. Missing DP values do not pass. Non-integer DP values are
rejected. Records may carry any ``FILTER`` value; the file is assumed to have
already been filtered as desired. Duplicate ``CHROM:POS`` records are combined
with OR semantics for each sample, but duplicate records must be contiguous, as
in a coordinate-sorted VCF.

Target BED
==========

``--targets`` accepts a BED file with at least three columns. Blank lines,
comments, and a ``chrom start end`` header are ignored. Coordinates must be
non-negative, 0-based, half-open intervals where ``start < end``.

In BAM/CRAM mode, targets are sorted and merged with ``bedtools`` before use.
In VCF mode, overlapping target intervals are merged internally before bases are
tested for containment.
