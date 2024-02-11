#!/bin/bash

set -ex

for cell_line in DMS53 HCC1187 HCC1395 K562 KIJK MJ RT112 SKBR3 VCAP; do
    ../scripts/examine_fusion_brkpt_3prime_read_dist.py \
        --fusion_predictions ${cell_line}.ctat-LR-fusion.fusion_predictions.tsv.gz \
        --gff3_read_alignments ${cell_line}.LR-FI.mm2.gff3.gz \
        --output_prefix ${cell_line}.3prime_brkpt_read_lengths
    echo done $cell_line
done

