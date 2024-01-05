#!/bin/bash


set -ex

####################################################################
# be sure to have first run:
#   make incl_extra_starF_supported_in_truth_by_Illum
# so results are based on including the Illumina-supported entries.
####################################################################


scp login://home/unix/bhaas/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/cancer_cell_lines/preds.collected.gencode_mapped.wAnnot.filt .

scp login://home/unix/bhaas/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/cancer_cell_lines/__min_2_agree/min_2.okPara_ignoreUnsure.results.scored .

scp login://home/unix/bhaas/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/cancer_cell_lines/__min_2_agree/min_2.okPara_ignoreUnsure.results.scored.ROC .

scp login://home/unix/bhaas/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/cancer_cell_lines/__min_2_agree/preds.collected.gencode_mapped.wAnnot.filt.proxy_assignments.byProgAgree.min_2.truth_set .

scp login://home/unix/bhaas/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/cancer_cell_lines/__min_2_agree/preds.collected.gencode_mapped.wAnnot.filt.proxy_assignments.byProgAgree.min_2.unique_set .

scp login://home/unix/bhaas/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/cancer_cell_lines/__min_2_agree/preds.collected.gencode_mapped.wAnnot.filt.proxy_assignments.byProgAgree.min_2.unsure_set .
