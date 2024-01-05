#!/bin/bash

set -ex

scp login://home/unix/bhaas/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/SeraCareFusions/__analyze_allow_reverse/combined_results.ROC.tsv seracarefusion.allow_rev.combined_results.ROC.tsv

scp login://home/unix/bhaas/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/SeraCareFusions/__analyze_allow_reverse/ISO-seq/fusion_preds.txt.scored Iso-seq.fusion_preds.txt.scored

scp login://home/unix/bhaas/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/SeraCareFusions/__analyze_allow_reverse/MAS-seq-L1/fusion_preds.txt.scored MAS-seq-L1.fusion_preds.txt.scored

scp login://home/unix/bhaas/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/SeraCareFusions/__analyze_allow_reverse/MAS-seq-L2/fusion_preds.txt.scored MAS-seq-L2.fusion_preds.txt.scored

gzip -f *scored

