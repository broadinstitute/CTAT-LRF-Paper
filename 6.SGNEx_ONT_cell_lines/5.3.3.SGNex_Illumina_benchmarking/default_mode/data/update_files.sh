#!/bin/bash

set -ex


rm -f *.gz

cp ~/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/SGNex_ONT/__valid_set_TP_uniq_FP/*set .

cp ~/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/SGNex_ONT/__valid_set_TP_uniq_FP/*ignoreUnsure.results.scored.ROC .

cp ~/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/SGNex_ONT/__valid_set_TP_uniq_FP/*ignoreUnsure.results.scored.PR* .

cp ~/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/SGNex_ONT/__valid_set_TP_uniq_FP/*ignoreUnsure.results.scored .

gzip *scored

