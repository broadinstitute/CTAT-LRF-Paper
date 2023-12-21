version 1.0


workflow LongGF {
    input {
        File fastq
        String sample_id
        String platform
        File fasta
        File gtf
    }

    call LongGFTask {
        input:
            input_fastq=fastq,
            sample_id=sample_id,
            platform=platform,
            fasta=fasta,
            gtf=gtf,
    }

    output {
        File confident_fusion = LongGFTask.output_confident_fusion
    }
}


task LongGFTask {
    input {
        File input_fastq
        String sample_id
        String platform
        File fasta
        File gtf

        String docker_image="trinityctat/longgf"
        String assembly="GRCh38"
        Int preemptible=3
        Int boot_disk_size=10
        Int disk_space=20
        Int cpu = 4
        Int mem = 32
    }

    command <<<
        minimap2 -t ~{cpu} -ax splice ~{fasta} ~{input_fastq} | samtools sort -n -o ~{sample_id}_sort.bam
        LongGF ~{sample_id}_sort.bam ~{gtf} 100 50 100 0 0 1 0 >~{sample_id}_longgf.out
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
        File output_confident_fusion = "~{sample_id}_longgf.out"
        File bam = "~{sample_id}_sort.bam"
    }
}
