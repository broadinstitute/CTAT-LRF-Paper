version 1.0


workflow JAFFAL {
    input {
        File fastq
        String sample_id
    }

    call JAFFALTask {
        input:
            input_fastq=fastq,
            sample_id=sample_id,
    }

    output {
        File confident_fusion = JAFFALTask.output_confident_fusion
        File confident_fusion_fa = JAFFALTask.output_confident_fusion_fa
        File output_barcodes_txt = JAFFALTask.output_barcodes_txt
        File output_transcriptome_paf = JAFFALTask.output_transcriptome_paf
        File output_genome_paf = JAFFALTask.output_genome_paf
        File output_genome_psl = JAFFALTask.output_genome_psl
        File output_fusion_fasta = JAFFALTask.output_fusion_fasta
        File output_reads = JAFFALTask.output_reads
        File output_summary = JAFFALTask.output_summary
        File output_log = JAFFALTask.output_log
    }
}


task JAFFALTask {
    input {
        File input_fastq
        String sample_id

        String docker_image="trinityctat/jaffal"
        String assembly="GRCh38"
        Int preemptible=3
        Int boot_disk_size=20
        Int disk_space=30
        Int cpu = 8
        Int mem = 32
    }

    command <<<

        /opt/JAFFA/tools/bin/bpipe run -n ~{cpu} /opt/JAFFA/JAFFAL.groovy ~{input_fastq}
        du -sh jaffa_results.csv jaffa_results.fasta
        echo $(pwd)

        mv jaffa_results.csv ~{sample_id}_jaffal_results.csv
        mv jaffa_results.fasta ~{sample_id}_jaffal_results.fasta

        label=$(basename ~{input_fastq})
        label=${label/.gz/}
        ls ${label}

        # from transcriptome alignment
        cp -r ${label}/${label}.paf ~{sample_id}.paf
        cp -r ${label}/${label}.txt ~{sample_id}_jaffal_results.txt
        cp -r ${label}/${label}.reads ~{sample_id}.reads

        cp -r ${label}/${label}.fusions.fa ~{sample_id}_fusions.fa
        # from genome alignment of fusion fasta
        cp -r ${label}/${label}_genome.paf ~{sample_id}_genome.paf
        cp -r ${label}/${label}_genome.psl ~{sample_id}_genome.psl

        # input for final output fusion table
        cp -r ${label}/${label}.summary ~{sample_id}.summary
        cp -r commandlog.txt ~{sample_id}_commandlog.txt

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
        File output_confident_fusion = "~{sample_id}_jaffal_results.csv"
        File output_confident_fusion_fa = "~{sample_id}_jaffal_results.fasta"
        File output_barcodes_txt = "~{sample_id}_jaffal_results.txt"
        File output_transcriptome_paf = "~{sample_id}.paf"
        File output_genome_paf = "~{sample_id}_genome.paf"
        File output_genome_psl = "~{sample_id}_genome.psl"
        File output_fusion_fasta = "~{sample_id}_fusions.fa"
        File output_reads = "~{sample_id}.reads"
        File output_summary = "~{sample_id}.summary"
        File output_log = "~{sample_id}_commandlog.txt"
    }
}
