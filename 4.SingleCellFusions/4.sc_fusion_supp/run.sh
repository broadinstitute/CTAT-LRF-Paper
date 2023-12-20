#!/bin/bash -ex

# first install https://github.com/ans4013/ScisorWiz
# mamba create -n r -c conda-forge r-devtools
# devtools::install_github('ans4013/ScisorWiz',build_vignettes = FALSE)
# then replace the error-prone script 
# /home/ubuntu/mambaforge/envs/r/lib/R/library/ScisorWiz/python/ClusterByIsoform_gffInput.py
# with ClusterByIsoform_gffInput.py in this folder

# short read
python parse_bam_readinfo_gff.py -s short -r igv_M132TS.NUTM2A-AS1--RP11-203L2.4/M132TS.FI.Illum10x.NUTM2A-AS1--RP11-203L2.4.bam.gtf -g igv_M132TS.NUTM2A-AS1--RP11-203L2.4/igv.annot.gtf > readinfo.gff
gzip readinfo.gff
gzip readstogenes_short
Rscript ScisorWiz_short.R

# long read
#perl igv_M132TS.NUTM2A-AS1--RP11-203L2.4/CTAT-LR-fusion/util/SAM_to_gxf.pl --sam igv_M132TS.NUTM2A-AS1--RP11-203L2.4/M132TS.NUTM2A-AS1--RP11-203L2.4.bam --format gtf > igv_M132TS.NUTM2A-AS1--RP11-203L2.4/M132TS.NUTM2A-AS1--RP11-203L2.4.bam.gtf
python parse_bam_readinfo_gff.py -r igv_M132TS.NUTM2A-AS1--RP11-203L2.4/M132TS.NUTM2A-AS1--RP11-203L2.4.bam.gtf -g igv_M132TS.NUTM2A-AS1--RP11-203L2.4/igv.annot.gtf > readinfo_long.gff
gzip readinfo_long.gff
gzip readstogenes_long
Rscript ScisorWiz_long.R
