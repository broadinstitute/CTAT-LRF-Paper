
mkdir arriba_output

singularity exec -B `pwd`/arriba_output:/output -B `pwd`/MDL_163_SeraSeqFusion_1.bam.left.fastq.gz:/read1.fastq.gz -B `pwd`/MDL_163_SeraSeqFusion_1.bam.right.fastq.gz:/read2.fastq.gz  -B /broad/hptmp/bhaas/ArribaReferences:/references ../../arriba-2.4.0.simg arriba.sh
