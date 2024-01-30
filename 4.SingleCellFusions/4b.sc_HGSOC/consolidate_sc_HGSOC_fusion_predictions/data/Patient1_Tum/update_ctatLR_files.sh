#!/bin/bash

set -ex

PREFIX="Patient1_Tum"
CTAT_LR_DIR="/seq/RNASEQ/FUSION_LONGREADS/Dondi/ORGANIZED_OUTPUTS/Patient1_Tum"


# files to be updated:

# -rw-r--r--  1 bhaas  staff   70943 Jan  2 11:22 Patient1_Om.ctat-LR-fusion.fusion_predictions.abridged.tsv.gz
# -rw-r--r--  1 bhaas  staff  161418 Jan  2 11:22 Patient1_Om.ctat-LR-fusion.fusion_predictions.tsv.gz

# -rw-r--r--  1 bhaas  staff   73807 Jan  2 11:22 Patient1_Om.ctatLRF_FI.fusion_predictions.abridged.tsv.gz
# -rw-r--r--  1 bhaas  staff  522853 Jan  2 11:22 Patient1_Om.ctatLRF_FI.fusion_predictions.tsv.gz




# CTAT-LR (no FI)

scp login:/${CTAT_LR_DIR}/PacBio/ctat_LR_fusion_outdir/ctat-LR-fusion.fusion_predictions.abridged.tsv ${PREFIX}.ctat-LR-fusion.fusion_predictions.abridged.tsv

scp login:/${CTAT_LR_DIR}/PacBio/ctat_LR_fusion_outdir/ctat-LR-fusion.fusion_predictions.tsv ${PREFIX}.ctat-LR-fusion.fusion_predictions.tsv


# CTAT-LRF_FI

scp login:/${CTAT_LR_DIR}/ctatLRF_pb_and_illum_together/ctat_LR_fusion_outdir/ctat-LR-fusion.fusion_predictions.abridged.tsv ${PREFIX}.ctatLRF_FI.fusion_predictions.abridged.tsv

scp login:/${CTAT_LR_DIR}/ctatLRF_pb_and_illum_together/ctat_LR_fusion_outdir/ctat-LR-fusion.fusion_predictions.tsv ${PREFIX}.ctatLRF_FI.fusion_predictions.tsv

gzip -f *tsv



