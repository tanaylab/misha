#!/bin/sh
# Rebuilds the BAM fixtures of the gtrack.import_mappedseq tests. Needs samtools and python3.
# The same files are in misha (tests/testthat/fixtures/bam) and pymisha (tests/fixtures/bam).
#
#   <name>.bam          samtools view -b of <name>.src.sam (python3 gen_sam.py . writes them)
#   split.bam           single.bam re-compressed into small BGZF blocks, with empty blocks among
#                       them, so that records and the header span blocks
#   edge.bam            records samtools cannot write from SAM text (a mapped record without a
#                       CIGAR, CIGAR operation codes SAM has no letter for, ...)
#   bad-*.bam           files htslib refuses
#   <name>.view.sam.gz  samtools view -h of <name>.bam: the text misha imported through samtools
set -e
cd "$(dirname "$0")"
python3 gen_sam.py .
for f in default single paired empty; do
    samtools view --no-PG -b -o $f.bam $f.src.sam
done
python3 bamtools.py reblock single.bam split.bam 7
python3 bamtools.py edge edge.bam
python3 bamtools.py corrupt .
for f in default single paired empty edge; do
    samtools view --no-PG -h $f.bam | gzip -9 -n > $f.view.sam.gz
done
rm -f single.src.sam paired.src.sam empty.src.sam
