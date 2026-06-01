Output files
************

Final outputs
=============

Each successful run writes two files to ``--out``:

.. code-block:: text

   sprite.bed.gz
   sprite.bed.gz.tbi

The population count mask is bgzip-compressed and indexed with
``tabix -p bed``. Use ``--output-prefix`` to choose a different filename
stem; ``.bed.gz`` is always appended.

BED columns
===========

The mask starts with two comment-prefixed header lines:

.. code-block:: text

   #sprite_mask_metadata  {"columns":["chrom","start","end","GBR","YRI"],...}
   #chrom  start  end  GBR  YRI

Data rows contain:

**chrom**
    Chromosome or contig name.

**start**
    0-based BED start coordinate.

**end**
    0-based BED end coordinate, exclusive.

**population columns**
    Number of samples in each population with depth greater than or equal to
    ``--min-dp`` over the interval.

Adjacent bases with identical population counts are collapsed. Intervals where
all counts are zero are omitted.

Metadata
========

The ``#sprite_mask_metadata`` line is JSON. It includes:

* ``sprite_mask_version``
* ``format``
* ``columns``
* ``coordinate_system``
* ``zero_count_intervals_omitted``
* ``min_dp``
* ``sample_count``
* ``populations``
* ``population_columns``
* ``population_sample_counts``
* input paths such as ``samples_path``, ``popfile``, ``all_sites_vcf``, and
  ``mask_bed`` when applicable
* ``variants_vcf``, ``threshold_sources``, and
  ``variant_vcf_threshold_estimates`` when ``--variants-vcf`` is used in
  BAM/CRAM mode

Sparse interpretation
=====================

The mask is sparse by design. A missing interval means zero passing samples
in every population — not that the interval was skipped or that counts are
unknown.

When ``--variants-vcf`` excludes an indel, structural variant, breakend, or
multi-nucleotide polymorphism span in BAM/CRAM mode, that span is removed from
all sample pass BEDs. In the final sparse BED, this is represented the same
way as any other all-zero region: no data row is written for the span.

Intermediate files
==================

When ``--keep-work`` is set, BAM/CRAM mode retains files like:

.. code-block:: text

   mask.3col.bed
   mask.3col.sorted.merged.bed
   <sample>.d<min-dp>.quantized.bed.gz
   <sample>.d<min-dp>.quantized.bed.gz.csi
   <sample>.d<min-dp>.mosdepth.summary.txt
   <sample>.d<min-dp>.mosdepth.global.dist.txt
   <sample>.d<min-dp>.mosdepth.stderr.log
   <sample>.d<min-dp>.pass.bed
   <sample>.d<min-dp>.pass.targets.bed
   variants_vcf.excluded.raw.bed
   variants_vcf.excluded.sorted.merged.bed
   <sample>.d<min-dp>.pass.variants.bed
   <sample>.d<min-dp>.pass.targets.variants.bed
   cohort.d<min-dp>.multiinter.tsv
   cohort.d<min-dp>.population_count_quantized.bed

The ``variants_vcf.*`` and ``*.variants.bed`` files are present only when
``--variants-vcf`` finds non-SNP exclusion intervals.

VCF mode retains the uncompressed population count mask in the work directory.
