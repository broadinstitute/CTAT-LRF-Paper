#!/bin/bash

set -ex

PREFIX="Patient3_Tum"
CTAT_LR_DIR="/seq/RNASEQ/FUSION_LONGREADS/CTAT_LR_Fusion-paper/tumor_single_cells/HGSOC/HGSOC.Patient3"

# files to be updated:

# -rw-r--r--  1 bhaas  staff   70943 Jan  2 11:22 Patient1_Om.ctat-LR-fusion.fusion_predictions.abridged.tsv.gz
# -rw-r--r--  1 bhaas  staff  161418 Jan  2 11:22 Patient1_Om.ctat-LR-fusion.fusion_predictions.tsv.gz

# -rw-r--r--  1 bhaas  staff   73807 Jan  2 11:22 Patient1_Om.ctatLRF_FI.fusion_predictions.abridged.tsv.gz
# -rw-r--r--  1 bhaas  staff  522853 Jan  2 11:22 Patient1_Om.ctatLRF_FI.fusion_predictions.tsv.gz




# CTAT-LR (no FI)

scp login:/${CTAT_LR_DIR}/HGSOC.Patient3.Tum/HGSOC.Patient3.Tum.v0.12.0.LRF/ctat-LR-fusion.fusion_predictions.abridged.tsv.gz ${PREFIX}.ctat-LR-fusion.fusion_predictions.abridged.tsv.gz

scp login:/${CTAT_LR_DIR}/HGSOC.Patient3.Tum/HGSOC.Patient3.Tum.v0.12.0.LRF/ctat-LR-fusion.fusion_predictions.tsv.gz ${PREFIX}.ctat-LR-fusion.fusion_predictions.tsv.gz

# CTAT-LRF_FI

scp login:/${CTAT_LR_DIR}/HGSOC.Patient3.Tum/HGSOC.Patient3.Tum.v0.12.0.LRF_FI/ctat-LR-fusion.fusion_predictions.abridged.tsv.gz ${PREFIX}.ctatLRF_FI.fusion_predictions.abridged.tsv.gz

scp login:/${CTAT_LR_DIR}/HGSOC.Patient3.Tum/HGSOC.Patient3.Tum.v0.12.0.LRF_FI/ctat-LR-fusion.fusion_predictions.tsv.gz ${PREFIX}.ctatLRF_FI.fusion_predictions.tsv.gz

