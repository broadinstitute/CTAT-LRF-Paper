#!/bin/bash

set -ex

scp login://home/unix/bhaas/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/cancer_cell_lines/preds.collected.gencode_mapped.wAnnot .

gzip -f preds.collected.gencode_mapped.wAnnot

scp login://home/unix/bhaas/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/cancer_cell_lines/preds.collected.gencode_mapped.wAnnot.filt.pass .



## min 2 agree analysis

X=2

scp login://home/unix/bhaas/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/cancer_cell_lines/__min_${X}_agree/min_${X}.okPara_ignoreUnsure.results.scored .

scp login://home/unix/bhaas/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/cancer_cell_lines/__min_${X}_agree/min_${X}.okPara_ignoreUnsure.results.scored.ROC .

scp login://home/unix/bhaas/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/cancer_cell_lines/__min_${X}_agree/preds.collected.gencode_mapped.wAnnot.filt.pass.proxy_assignments.byProgAgree.min_${X}.truth_set .

scp login://home/unix/bhaas/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/cancer_cell_lines/__min_${X}_agree/preds.collected.gencode_mapped.wAnnot.filt.pass.proxy_assignments.byProgAgree.min_${X}.unique_set .

scp login://home/unix/bhaas/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/cancer_cell_lines/__min_${X}_agree/preds.collected.gencode_mapped.wAnnot.filt.pass.proxy_assignments.byProgAgree.min_${X}.unsure_set .


scp login://home/unix/bhaas/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/cancer_cell_lines/__min_${X}_agree/min_${X}.ignoreUnsure.results.scored .

scp login://home/unix/bhaas/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/cancer_cell_lines/__min_${X}_agree/min_${X}.ignoreUnsure.results.scored.ROC .


## min 3 agree analysis

X=3

scp login://home/unix/bhaas/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/cancer_cell_lines/__min_${X}_agree/min_${X}.okPara_ignoreUnsure.results.scored .

scp login://home/unix/bhaas/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/cancer_cell_lines/__min_${X}_agree/min_${X}.okPara_ignoreUnsure.results.scored.ROC .

scp login://home/unix/bhaas/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/cancer_cell_lines/__min_${X}_agree/preds.collected.gencode_mapped.wAnnot.filt.pass.proxy_assignments.byProgAgree.min_${X}.truth_set .

scp login://home/unix/bhaas/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/cancer_cell_lines/__min_${X}_agree/preds.collected.gencode_mapped.wAnnot.filt.pass.proxy_assignments.byProgAgree.min_${X}.unique_set .

scp login://home/unix/bhaas/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/cancer_cell_lines/__min_${X}_agree/preds.collected.gencode_mapped.wAnnot.filt.pass.proxy_assignments.byProgAgree.min_${X}.unsure_set .


scp login://home/unix/bhaas/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/cancer_cell_lines/__min_${X}_agree/min_${X}.ignoreUnsure.results.scored .

scp login://home/unix/bhaas/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/cancer_cell_lines/__min_${X}_agree/min_${X}.ignoreUnsure.results.scored.ROC .
