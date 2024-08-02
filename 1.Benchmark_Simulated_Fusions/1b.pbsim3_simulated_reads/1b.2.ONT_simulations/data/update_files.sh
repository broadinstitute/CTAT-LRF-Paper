#!/bin/bash

set -ex

# breakpoint_maxF1_data.tsv max_F1_summary.tsv


cp ~/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/simulated_data/our_pbsim3_sims/ONT_pbsim3_part5/max_F1_summary.tsv .

cp ~/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/simulated_data/our_pbsim3_sims/ONT_pbsim3_part5/breakpoint_maxF1_data.tsv .


cp ~/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/simulated_data/our_pbsim3_sims/ONT_pbsim3_part5/all_combined_results.PR_AUC.tsv .

cp ~/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/simulated_data/our_pbsim3_sims/ONT_pbsim3_part5/mean_PR_AUC_summary.tsv .

cp ~/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/simulated_data/our_pbsim3_sims/ONT_pbsim3_part5/breakpoint_all_PR_AUC_results.tsv .

cp ~/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/simulated_data/our_pbsim3_sims/ONT_pbsim3_part5/breakpoint_all_mean_PR_AUC_data.tsv .



#scp login://home/unix/bhaas/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/simulated_data/ONT_pbsim3_part5/max_F1_summary.tsv .

#scp login://home/unix/bhaas/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/simulated_data/ONT_pbsim3_part5/breakpoint_maxF1_data.tsv .

