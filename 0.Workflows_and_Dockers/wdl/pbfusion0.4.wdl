version 1.0


workflow pbfusion {
    input {
        File fastq
        String sample_id
        File mmi
        File gtf_bin
        String platform
        String docker_image="trinityctat/pbfusion:v0.4.0"
        String version="v0.4.0"
        String assembly="GRCh38"
        Int preemptible=3
        Int boot_disk_size=10
        Int disk_space=20
        Int cpu = 10
        Int mem = 32
        Int threshold=0
        # only apply to 0.3.0strand version
        String strand_option=" --strand-insensitive "
    }

    call pbfusionTask {
        input:
            input_fastq=fastq,
            sample_id=sample_id,
            mmi=mmi,
            gtf_bin=gtf_bin,
            platform=platform,
            docker_image=docker_image,
            assembly=assembly,
            preemptible=preemptible,
            boot_disk_size=boot_disk_size,
            disk_space=disk_space,
            cpu=cpu,
            mem=mem,
            threshold=threshold,
            strand_option=strand_option,
            version=version
    }

    output {
        File confident_fusion = pbfusionTask.output_confident_fusion
        File bam = pbfusionTask.bam
    }
}


task pbfusionTask {
    input {
        File input_fastq
        String sample_id
        File mmi
        File gtf_bin
        String platform

        String docker_image="trinityctat/pbfusion:v0.4.0"
        String version="v0.4.0"
        String assembly="GRCh38"
        Int preemptible=3
        Int boot_disk_size=10
        Int disk_space=20
        Int cpu = 10
        Int mem = 32
        Int threshold=0
        String strand_option=" --strand-insensitive "
    }

    command <<<
        if [[ ~{platform} == 'pacbio' ]]; then
            seqtk seq -A ~{input_fastq}  | awk '{if (/^>/){print $1"-molecule"} else {print}}' > ~{sample_id}.fa
            pbmm2 align --num-threads ~{cpu} --preset ISOSEQ --sort ~{mmi} ~{sample_id}.fa ~{sample_id}.bam
        else
            seqtk seq -A ~{input_fastq}  | awk '{if (/^>/){print $1"-ccs"} else {print}}' > ~{sample_id}.fa
            pbmm2 align --num-threads ~{cpu} --preset CCS --sort ~{mmi} ~{sample_id}.fa ~{sample_id}.bam
        fi

        if [[ ~{version} == 'v0.1.0' ]]; then
            pbfusion --reference-gtf ~{gtf_bin} --bam ~{sample_id}.bam --output-prefix ~{sample_id}_pbfusion --threads ~{cpu} -v --min-coverage ~{threshold} 
        elif [[ ~{version} == 'v0.4.0' ]]; then
            pbfusion discover --gtf ~{gtf_bin} --output-prefix ~{sample_id}_pbfusion --threads ~{cpu} -v --min-coverage ~{threshold} ~{sample_id}.bam
        else
            pbfusion discover --gtf ~{gtf_bin} --bam ~{sample_id}.bam --output-prefix ~{sample_id}_pbfusion --threads ~{cpu} -v --min-coverage ~{threshold} ~{strand_option}
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
        File output_confident_fusion = "~{sample_id}_pbfusion.breakpoints.groups.bed"
        File bam = "~{sample_id}.bam"
    }
}
