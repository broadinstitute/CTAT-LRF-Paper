#!/bin/bash

set -ex

scp login://home/unix/bhaas/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/cancer_cell_lines/preds.collected.gencode_mapped.wAnnot .

gzip -f preds.collected.gencode_mapped.wAnnot

scp login://home/unix/bhaas/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/cancer_cell_lines/preds.collected.gencode_mapped.wAnnot.filt.pass .

scp login://home/unix/bhaas/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/cancer_cell_lines/__min_2_agree/min_2.okPara_ignoreUnsure.results.scored .

scp login://home/unix/bhaas/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/cancer_cell_lines/__min_2_agree/min_2.okPara_ignoreUnsure.results.scored.ROC .

scp login://home/unix/bhaas/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/cancer_cell_lines/__min_2_agree/preds.collected.gencode_mapped.wAnnot.filt.pass.proxy_assignments.byProgAgree.min_2.truth_set .

scp login://home/unix/bhaas/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/cancer_cell_lines/__min_2_agree/preds.collected.gencode_mapped.wAnnot.filt.pass.proxy_assignments.byProgAgree.min_2.unique_set .

scp login://home/unix/bhaas/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/cancer_cell_lines/__min_2_agree/preds.collected.gencode_mapped.wAnnot.filt.pass.proxy_assignments.byProgAgree.min_2.unsure_set .


scp login://home/unix/bhaas/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/cancer_cell_lines/__min_2_agree/min_2.ignoreUnsure.results.scored .

scp login://home/unix/bhaas/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/cancer_cell_lines/__min_2_agree/min_2.ignoreUnsure.results.scored.ROC .


