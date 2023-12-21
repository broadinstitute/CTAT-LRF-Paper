version 1.0


workflow FusionSeeker {
    input {
        File fastq
        String sample_id
        String platform
        File fasta
        File gtf
        String docker_image="qianqin/fusion_seeker:latest"
        String assembly="GRCh38"
        Int preemptible=3
        Int boot_disk_size=10
        Int disk_space=20
        Int cpu = 4
        Int mem = 32
        Int readthresh = 1
    }

    call FusionSeekerTask {
        input:
            input_fastq=fastq,
            sample_id=sample_id,
            platform=platform,
            fasta=fasta,
            gtf=gtf,
            docker_image=docker_image,
            assembly=assembly,
            preemptible=preemptible,
            boot_disk_size=boot_disk_size,
            disk_space=disk_space,
            cpu=cpu,
            mem=mem,
            readthresh=readthresh
    }

    output {
        File confident_fusion = FusionSeekerTask.output_confident_fusion
    }
}


task FusionSeekerTask {
    input {
        File input_fastq
        String sample_id
        String platform
        File fasta
        File gtf

        String docker_image="qianqin/fusion_seeker:latest"
        String assembly="GRCh38"
        Int preemptible=3
        Int boot_disk_size=10
        Int disk_space=20
        Int cpu = 4
        Int mem = 32
        Int readthresh = 1
    }

    command <<<
        if [[ ~{platform} == 'pacbio' ]]; then
            minimap2 -t ~{cpu} -ax splice:hq ~{fasta} ~{input_fastq} | samtools sort -o ~{sample_id}.bam
            samtools index ~{sample_id}.bam
            fusionseeker --bam ~{sample_id}.bam -o ~{sample_id}_output --datatype isoseq --ref ~{fasta} --gtf ~{gtf} -s ~{readthresh}
        else
            minimap2 -t ~{cpu} -ax splice ~{fasta} ~{input_fastq} | samtools sort -o ~{sample_id}.bam
            samtools index ~{sample_id}.bam 
            fusionseeker --bam ~{sample_id}.bam -o ~{sample_id}_output --datatype nanopore --ref ~{fasta} --gtf ~{gtf} -s ~{readthresh}
        fi
    >>>

    runtime {
        disks: "local-disk ~{disk_space} HDD"
        memory: "~{mem} GB"
        cpu: cpu
        preemptible: preemptible
        bootDiskSizeGb: boot_disk_size
        docker: docker_image
    }

    output {
        File output_confident_fusion = "~{sample_id}_output/confident_genefusion.txt"
        File bam = "~{sample_id}.bam"
    }
}
