#!/bin/bash

set -ex


rm -f *.gz

cp ~/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/DepMap_Cell_Lines/preds.collected.gencode_mapped.wAnnot.filt.pass.proxy_assignments .

cp  ~/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/DepMap_Cell_Lines/preds.collected.gencode_mapped.wAnnot.filt.pass.proxy_assignments.byProgAgree .

gzip -f preds.collected.gencode_mapped*



DIRNAME="__illum_TP_uniq_FP.starF"


cp ~/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/DepMap_Cell_Lines/${DIRNAME}/preds.collected.gencode_mapped.wAnnot.filt.pass.proxy_assignments.byProgAgree.illum_agree.truth_set .


cp ~/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/DepMap_Cell_Lines/${DIRNAME}/preds.collected.gencode_mapped.wAnnot.filt.pass.proxy_assignments.byProgAgree.nonunique.unsure_set .


cp ~/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/DepMap_Cell_Lines/${DIRNAME}/eval_illum_supported.ignoreUnsure.results.scored .

cp ~/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/DepMap_Cell_Lines/${DIRNAME}/eval_illum_supported.ignoreUnsure.results.scored.ROC .

cp ~/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/DepMap_Cell_Lines/${DIRNAME}/eval_illum_supported.ignoreUnsure.results.scored.PR .

cp ~/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/DepMap_Cell_Lines/${DIRNAME}/eval_illum_supported.ignoreUnsure.results.scored.PR.AUC .



cp ~/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/DepMap_Cell_Lines/${DIRNAME}/eval_illum_supported.okPara_ignoreUnsure.results.scored .

cp ~/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/DepMap_Cell_Lines/${DIRNAME}/eval_illum_supported.okPara_ignoreUnsure.results.scored.ROC .

cp ~/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/DepMap_Cell_Lines/${DIRNAME}/eval_illum_supported.okPara_ignoreUnsure.results.scored.PR .

cp ~/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/DepMap_Cell_Lines/${DIRNAME}/eval_illum_supported.okPara_ignoreUnsure.results.scored.PR.AUC .












