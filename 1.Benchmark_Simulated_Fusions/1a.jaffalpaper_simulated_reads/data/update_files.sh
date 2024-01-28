#!/bin/bash

set -ex


# allow_reverse.combined_results.ROC.tsv strict.combined_results.ROC.tsv

scp login://home/unix/bhaas/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/simulated_data/sim_jaffal/__analyze_strict/combined_results.ROC.tsv strict.combined_results.ROC.tsv

scp login://home/unix/bhaas/GITHUB/CTAT_FUSIONS/LR-FusionBenchmarking/simulated_data/sim_jaffal/__analyze_allow_reverse/combined_results.ROC.tsv allow_reverse.combined_results.ROC.tsv

