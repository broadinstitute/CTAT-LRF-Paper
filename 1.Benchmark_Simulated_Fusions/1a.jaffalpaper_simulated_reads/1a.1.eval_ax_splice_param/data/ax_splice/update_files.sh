#!/bin/bash

set -ex


# allow_reverse.combined_results.ROC.tsv strict.combined_results.ROC.tsv


cp ~/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/simulated_data/sim_jaffal/__analyze_strict_and_paralogs/combined_results.ROC.tsv strict_allow_paralogs.combined_results.ROC.tsv

cp ~/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/simulated_data/sim_jaffal/__analyze_allow_rev_and_paralogs/combined_results.ROC.tsv allow_rev_and_paralogs.combined_results.ROC.tsv


cp  ~/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/simulated_data/sim_jaffal/__analyze_strict_and_paralogs/combined_results.PR_AUC.tsv strict_allow_paralogs.combined_results.PR_AUC.tsv

cp ~/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/simulated_data/sim_jaffal/__analyze_allow_rev_and_paralogs/combined_results.PR_AUC.tsv allow_rev_and_paralogs.combined_results.PR_AUC.tsv   



#scp login://home/unix/bhaas/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/simulated_data/sim_jaffal/__analyze_strict/combined_results.ROC.tsv strict.combined_results.ROC.tsv

#scp login://home/unix/bhaas/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/simulated_data/sim_jaffal/__analyze_allow_reverse/combined_results.ROC.tsv allow_reverse.combined_results.ROC.tsv

