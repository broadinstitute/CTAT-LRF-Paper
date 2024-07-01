#!/bin/bash

set -ex

cp  ~/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/DepMap_Cell_Lines/preds.collected.gencode_mapped.wAnnot .

gzip -f preds.collected.gencode_mapped.wAnnot

cp  ~/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/DepMap_Cell_Lines/preds.collected.gencode_mapped.wAnnot.filt.pass .



## min 2 agree analysis

X=2

cp  ~/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/DepMap_Cell_Lines/__min_${X}_agree/min_${X}.okPara_ignoreUnsure.results.scored .

cp  ~/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/DepMap_Cell_Lines/__min_${X}_agree/min_${X}.okPara_ignoreUnsure.results.scored.ROC .

cp  ~/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/DepMap_Cell_Lines/__min_${X}_agree/preds.collected.gencode_mapped.wAnnot.filt.pass.proxy_assignments.byProgAgree.min_${X}.truth_set .

cp  ~/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/DepMap_Cell_Lines/__min_${X}_agree/preds.collected.gencode_mapped.wAnnot.filt.pass.proxy_assignments.byProgAgree.min_${X}.unique_set .

cp  ~/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/DepMap_Cell_Lines/__min_${X}_agree/preds.collected.gencode_mapped.wAnnot.filt.pass.proxy_assignments.byProgAgree.min_${X}.unsure_set .


cp  ~/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/DepMap_Cell_Lines/__min_${X}_agree/min_${X}.ignoreUnsure.results.scored .

cp  ~/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/DepMap_Cell_Lines/__min_${X}_agree/min_${X}.ignoreUnsure.results.scored.ROC .


## min 3 agree analysis

X=3

cp  ~/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/DepMap_Cell_Lines/__min_${X}_agree/min_${X}.okPara_ignoreUnsure.results.scored .

cp  ~/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/DepMap_Cell_Lines/__min_${X}_agree/min_${X}.okPara_ignoreUnsure.results.scored.ROC .

cp  ~/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/DepMap_Cell_Lines/__min_${X}_agree/preds.collected.gencode_mapped.wAnnot.filt.pass.proxy_assignments.byProgAgree.min_${X}.truth_set .

cp  ~/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/DepMap_Cell_Lines/__min_${X}_agree/preds.collected.gencode_mapped.wAnnot.filt.pass.proxy_assignments.byProgAgree.min_${X}.unique_set .

cp  ~/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/DepMap_Cell_Lines/__min_${X}_agree/preds.collected.gencode_mapped.wAnnot.filt.pass.proxy_assignments.byProgAgree.min_${X}.unsure_set .


cp  ~/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/DepMap_Cell_Lines/__min_${X}_agree/min_${X}.ignoreUnsure.results.scored .

cp  ~/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/DepMap_Cell_Lines/__min_${X}_agree/min_${X}.ignoreUnsure.results.scored.ROC .
