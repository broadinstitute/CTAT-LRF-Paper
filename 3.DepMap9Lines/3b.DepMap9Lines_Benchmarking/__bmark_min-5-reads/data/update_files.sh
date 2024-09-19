#!/bin/bash

set -ex

## min X agree analysis


for X in 2 3 4; do

cp  ~/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/DepMap_Cell_Lines/__min_${X}_agree/preds.collected.gencode_mapped.wAnnot.filt.pass.proxy_assignments.byProgAgree.min_${X}.truth_set .

cp  ~/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/DepMap_Cell_Lines/__min_${X}_agree/preds.collected.gencode_mapped.wAnnot.filt.pass.proxy_assignments.byProgAgree.min_${X}.unique_set .

cp  ~/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/DepMap_Cell_Lines/__min_${X}_agree/preds.collected.gencode_mapped.wAnnot.filt.pass.proxy_assignments.byProgAgree.min_${X}.unsure_set .


cp  ~/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/DepMap_Cell_Lines/__min_${X}_agree/min_${X}.okPara_ignoreUnsure.results.scored .

cp  ~/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/DepMap_Cell_Lines/__min_${X}_agree/min_${X}.okPara_ignoreUnsure.results.scored.ROC .

cp ~/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/DepMap_Cell_Lines/__min_${X}_agree/min_${X}.okPara_ignoreUnsure.results.scored.PR .

cp ~/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/DepMap_Cell_Lines/__min_${X}_agree/min_${X}.okPara_ignoreUnsure.results.scored.PR.AUC .



cp  ~/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/DepMap_Cell_Lines/__min_${X}_agree/min_${X}.ignoreUnsure.results.scored .

cp  ~/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/DepMap_Cell_Lines/__min_${X}_agree/min_${X}.ignoreUnsure.results.scored.ROC .

cp  ~/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/DepMap_Cell_Lines/__min_${X}_agree/min_${X}.ignoreUnsure.results.scored.PR .

cp  ~/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/DepMap_Cell_Lines/__min_${X}_agree/min_${X}.ignoreUnsure.results.scored.PR.AUC .

done


