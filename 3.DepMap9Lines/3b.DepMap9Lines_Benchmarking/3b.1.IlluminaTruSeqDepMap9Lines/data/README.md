Run on Terra:
        https://app.terra.bio/#workspaces/methods-dev-lab/LR_Fusion_Benchmarking/job_history/c150ae25-fe49-4013-8421-c0ef023a73a8


example config:
    docker
"trinityctat/starfusion:1.10.1"
examine_coding_effect
false
extra_disk_space
10
fastq_disk_space_multiplier
3.25
fastq_pair_tar_gz
null
fusion_inspector
"validate"
genome_disk_space_multiplier
2.5
genome_plug_n_play_tar_gz
gs://ctat_genome_libs/__genome_libs_StarFv1.10/GRCh38_gencode_v22_CTAT_lib_Mar012021.plug-n-play.tar.gz
left_fq
gs://fc-d438f734-542a-445f-96e1-e07b4a7cddbc/DepMap/TruSeq/MDL_163_DMS53.bam.left.fastq.gz
memory
"50G"
num_cpu
12
preemptible
2
right_fq
gs://fc-d438f734-542a-445f-96e1-e07b4a7cddbc/DepMap/TruSeq/MDL_163_DMS53.bam.right.fastq.gz
sample_id
"DepMap_v1v2mrgd_DMS53"
use_ssd
true


