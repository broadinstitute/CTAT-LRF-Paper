version 1.0

struct RuntimeEnvironment {
    String docker_image
    Int preemptible
    Int boot_disk_size
    Int disk_space
    Int cpu
    Int mem 
}

workflow mergefusionWorkflow {
    meta {
        version: "v0.1"
        author: "Qian Qin"
    }

    input {
        File? fastq
        File? ubam
        String sample_id
        String readType
        File pbmmi
        File pbfusion_gtf_bin
        String platform
        String pbmm2_platform = 'pacbio'
        

        File genome_lib_tar_mm2_only
        File genome_lib_tar_with_STAR_idx
        Int min_per_id=90
        Int min_J = 1
        Int min_sumJS = 1    
        Int min_novel_junction_support = 1
        File? illumina_left_fq
        File? illumina_right_fq
      
        File fasta
        File gtf

        Float disk_space_multiplier = 3.0
        Int cpu = 10
        Int mem = 50
        Int disk_space
        Int preemptible
        Int boot_disk_size

        Boolean run_longgf = true
        Boolean run_fusionseeker = true
        Boolean run_jaffal = true
        Boolean run_pbfusion = true
        Boolean run_ctat = true

        File monitoring_script = "gs://broad-dsde-methods-monitoring/cromwell_monitoring_script.sh"
    }

    RuntimeEnvironment runtime_environment = {
        'docker_image': "trinityctat/longgf", 'preemptible': preemptible, 'boot_disk_size': boot_disk_size, 'disk_space': disk_space, 
        'cpu': cpu, 'mem': mem
    }

    RuntimeEnvironment runtime_environment2 = {
        'docker_image': "qianqin/fusion_seeker:latest", 'preemptible': preemptible, 'boot_disk_size': boot_disk_size, 'disk_space': disk_space, 
        'cpu': cpu, 'mem': mem
    }

    RuntimeEnvironment runtime_environment3 = {
        'docker_image': "trinityctat/jaffal", 'preemptible': preemptible, 'boot_disk_size': boot_disk_size, 'disk_space': disk_space, 
        'cpu': cpu, 'mem': mem
    }

    RuntimeEnvironment runtime_environment4 = {
        'docker_image': "trinityctat/pbfusion:v0.4.0", 'preemptible': preemptible, 'boot_disk_size': boot_disk_size, 'disk_space': disk_space, 
        'cpu': cpu, 'mem': mem
    }

    File genome_lib_tar = if defined(illumina_left_fq) then genome_lib_tar_with_STAR_idx else genome_lib_tar_mm2_only
    Int ctat_disk_space = ceil( (size(genome_lib_tar, "GB") + size(fastq, "GB") + size(ubam, "GB") + 2*size(illumina_left_fq, "GB") ) * disk_space_multiplier )
    RuntimeEnvironment runtime_environment5 = {
        'docker_image': "trinityctat/ctat_lr_fusion:latest", 'preemptible': preemptible, 'boot_disk_size': boot_disk_size, 'disk_space': ctat_disk_space, 
        'cpu': cpu, 'mem': mem
    }

    if ( defined(ubam) ) {
        call samtools_fastq {
            input:
                sample_id = sample_id,
                ubam = ubam,
                runtime_environment = runtime_environment,
                monitoring_script = monitoring_script
        }
    }

    File input_fastq = select_first([fastq, samtools_fastq.fastq])

    if ( !defined(input_fastq) ) {
        call raise_exception as error_input_data  { 
            input:
                msg = "No FASTQ input",
                runtime_environment = runtime_environment
        }
    }

    if ( defined(input_fastq) ) {
        if (run_ctat) {
            call CTAT_LR_FUSION_TASK as ctat_task {
               input:
                 sample_name=sample_id,
                 transcripts=input_fastq,
                 genome_lib_tar=genome_lib_tar,
                 min_per_id=min_per_id,
                 min_J=min_J,
                 min_sumJS=min_sumJS,    
                 min_novel_junction_support=min_novel_junction_support,
                 illumina_left_fq=illumina_left_fq,
                 illumina_right_fq=illumina_right_fq,

                 runtime_environment = runtime_environment5,
                 monitoring_script = monitoring_script
            }
        }
        if (run_pbfusion) {
            call seqtkTask {
                input:
                    input_fastq=input_fastq,
                    sample_id=sample_id,
                    monitoring_script = monitoring_script,
                    platform = platform,
                    runtime_environment = runtime_environment4
            }
            call pbmm2 {
                input:
                    input_fasta = seqtkTask.fa,
                    sample_id = sample_id,
                    pbmmi = pbmmi,
                    pbmm2_platform = pbmm2_platform,
                    monitoring_script = monitoring_script,
                    runtime_environment = runtime_environment4
            }
            call pbfusionTask {
                input:
                    bam = pbmm2.bam,
                    sample_id = sample_id,
                    pbfusion_gtf_bin = pbfusion_gtf_bin,
                    monitoring_script = monitoring_script,
                    runtime_environment = runtime_environment4
            }
        }

        if (run_jaffal) {
            call JAFFALTask {
                input:
                    input_fastq=input_fastq,
                    sample_id=sample_id,
                    monitoring_script = monitoring_script,
                    runtime_environment = runtime_environment3
            }
        }

        if ( run_longgf || run_fusionseeker ) {
            call minimap2Task {
                input:
                     fastq = input_fastq,
                     sample_id = sample_id,
                     gtf   = gtf,
                     fasta = fasta,
                     readType = readType,
                     runtime_environment = runtime_environment,
                     monitoring_script = monitoring_script
            }
        }
        if (run_longgf) {
            call LongGFTask {
                input:
                    bam = minimap2Task.bam,
                    sample_id = sample_id,
                    gtf = gtf,
                    runtime_environment = runtime_environment,
                    monitoring_script = monitoring_script
            }
        }

        if (run_fusionseeker) {
            call FusionSeekerTask {
                input:
                    bam = minimap2Task.bam,
                    sample_id = sample_id,
                    platform = platform,
                    gtf = gtf,
                    fasta = fasta,
                    runtime_environment = runtime_environment2,
                    monitoring_script = monitoring_script
            }
        }
    } 
    
    output {
        File? jaffal_confident_fusion = JAFFALTask.output_confident_fusion
        File? jaffal_confident_fusion_fa = JAFFALTask.output_confident_fusion_fa
        File? jaffal_output_barcodes_txt = JAFFALTask.output_barcodes_txt
        File? jaffal_output_transcriptome_paf = JAFFALTask.output_transcriptome_paf
        File? jaffal_output_genome_paf = JAFFALTask.output_genome_paf
        File? jaffal_output_genome_psl = JAFFALTask.output_genome_psl
        File? jaffal_output_fusion_fasta = JAFFALTask.output_fusion_fasta
        File? jaffal_output_reads = JAFFALTask.output_reads
        File? jaffal_output_summary = JAFFALTask.output_summary
        File? jaffal_output_log = JAFFALTask.output_log
        File? jaffal_log = JAFFALTask.monitoring_log

        File? ctat_fusion_report = ctat_task.fusion_report
        File? ctat_fusion_report_abridged = ctat_task.fusion_report_abridged

        File? ctat_prelim_fusion_report = ctat_task.prelim_fusion_report
        File? ctat_prelim_fusion_report_abridged = ctat_task.prelim_fusion_report_abridged

        File? ctat_fusion_report_html = ctat_task.fusion_report_html
        File? ctat_igv_tar = ctat_task.igv_tar
        File? ctat_monitoring_log = ctat_task.monitoring_log

        File? seqtk_log = seqtkTask.monitoring_log
        File? pbmm2_log = pbmm2.monitoring_log
        File? pbmm2_bam = pbmm2.bam
        File? pbfusion_log = pbfusionTask.monitoring_log
        File? pbfusion_fusion = pbfusionTask.output_confident_fusion

        File? longgf_confident_fusion = LongGFTask.output_confident_fusion
        File? fusionseeker_fusion = FusionSeekerTask.output_confident_fusion
        File? samtools_log = samtools_fastq.monitoring_log
        File? minimap2_log = minimap2Task.monitoring_log
        File? longgf_log = LongGFTask.monitoring_log
        File? fusionseeker_log = FusionSeekerTask.monitoring_log
    }
}


task samtools_fastq {
    input {
        String sample_id
        File? ubam
        RuntimeEnvironment runtime_environment

        File monitoring_script 
    }

    command <<< 
        bash ~{monitoring_script} > ~{sample_id}_monitoring.log &

        samtools fastq -@ ~{runtime_environment.cpu} ~{ubam} | gzip - > ~{sample_id}.fastq.gz
    >>>

    output {
        File? fastq = "~{sample_id}.fastq.gz"
        File? monitoring_log = "~{sample_id}_monitoring.log"
    }

    runtime {
        disks: "local-disk ~{runtime_environment.disk_space} HDD"
        memory: "~{runtime_environment.mem} GB"
        cpu: runtime_environment.cpu
        preemptible: runtime_environment.preemptible
        bootDiskSizeGb: runtime_environment.boot_disk_size
        docker: runtime_environment.docker_image
    }
}


task minimap2Task {
    input {
        File fastq
        String sample_id
        String readType
        File fasta
        File gtf

        RuntimeEnvironment runtime_environment
        File monitoring_script 
    }

    command <<<
        bash ~{monitoring_script} > ~{sample_id}_monitoring.log &

        # copy from https://github.com/broadinstitute/MDL-workflows/blob/main/LR-tools/minimap2_LR/minimap2_LR.wdl
        if [ "~{readType}" == "SplicedLongReads" ]; then
            minimap2_preset="splice"
        elif [ "~{readType}" == "ONTDirectRNA" ]; then
            minimap2_preset="splice -uf -k14"
        elif [ "~{readType}" == "PacBioIsoSeq" ]; then
            minimap2_preset="splice:hq "
        elif [ "~{readType}" == "PacBioIsoSeqUf" ]; then
            minimap2_preset="splice:hq -uf"
        elif [ "~{readType}" == "PacBioIsoSeqUb" ]; then
            minimap2_preset="splice:hq -ub"
        else
            echo "Invalid readType: ~{readType}"
            exit 1
        fi

        paftools.js gff2bed ~{gtf} > junc.bed
        minimap2 -t ~{runtime_environment.cpu} -ax ${minimap2_preset} --junc-bed junc.bed ~{fasta} ~{fastq} | samtools sort -n -o ~{sample_id}_sort.bam
    >>>

    runtime {
        disks: "local-disk ~{runtime_environment.disk_space} HDD"
        memory: "~{runtime_environment.mem} GB"
        cpu: runtime_environment.cpu
        preemptible: runtime_environment.preemptible
        bootDiskSizeGb: runtime_environment.boot_disk_size
        docker: runtime_environment.docker_image
    }

    output {
        File bam = "~{sample_id}_sort.bam"
        File monitoring_log = "~{sample_id}_monitoring.log"
    }
}

task LongGFTask {
    input {
        File? bam
        String sample_id
        File gtf
        RuntimeEnvironment runtime_environment
        File monitoring_script 
    }

    command <<<
        bash ~{monitoring_script} > ~{sample_id}_monitoring.log &
        LongGF ~{bam} ~{gtf} 100 50 100 0 0 1 0 >~{sample_id}_longgf.out
    >>>

    runtime {
        disks: "local-disk ~{runtime_environment.disk_space} HDD"
        memory: "~{runtime_environment.mem} GB"
        cpu: runtime_environment.cpu
        preemptible: runtime_environment.preemptible
        bootDiskSizeGb: runtime_environment.boot_disk_size
        docker: runtime_environment.docker_image
    }

    output {
        File output_confident_fusion = "~{sample_id}_longgf.out"
        File monitoring_log = "~{sample_id}_monitoring.log"
    }
}

task FusionSeekerTask {
    input {
        File? bam
        String sample_id
        String platform
        File gtf
        File fasta
        RuntimeEnvironment runtime_environment
        File monitoring_script 

        Int readthresh = 1
    }

    command <<<
        bash ~{monitoring_script} > ~{sample_id}_monitoring.log &
        samtools sort -o ~{sample_id}.bam ~{bam}
        samtools index ~{sample_id}.bam 
        if [[ ~{platform} == 'pacbio' ]]; then
            fusionseeker --bam ~{sample_id}.bam -o ~{sample_id}_output --datatype isoseq --ref ~{fasta} --gtf ~{gtf} -s ~{readthresh}
        else
            fusionseeker --bam ~{sample_id}.bam -o ~{sample_id}_output --datatype nanopore --ref ~{fasta} --gtf ~{gtf} -s ~{readthresh}
        fi
    >>>

    runtime {
        disks: "local-disk ~{runtime_environment.disk_space} HDD"
        memory: "~{runtime_environment.mem} GB"
        cpu: runtime_environment.cpu
        preemptible: runtime_environment.preemptible
        bootDiskSizeGb: runtime_environment.boot_disk_size
        docker: runtime_environment.docker_image
    }

    output {
        File? output_confident_fusion = "~{sample_id}_output/confident_genefusion.txt"
        File? monitoring_log = "~{sample_id}_monitoring.log"
    }
}

task JAFFALTask {
    input {
        File input_fastq
        String sample_id
        RuntimeEnvironment runtime_environment
        File monitoring_script 
    }

    command <<<
        bash ~{monitoring_script} > ~{sample_id}_monitoring.log &

        /opt/JAFFA/tools/bin/bpipe run -n ~{runtime_environment.cpu} /opt/JAFFA/JAFFAL.groovy ~{input_fastq}

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
        disks: "local-disk ~{runtime_environment.disk_space} HDD"
        memory: "~{runtime_environment.mem} GB"
        cpu: runtime_environment.cpu
        preemptible: runtime_environment.preemptible
        bootDiskSizeGb: runtime_environment.boot_disk_size
        docker: runtime_environment.docker_image
    }

    output {
        File? output_confident_fusion = "~{sample_id}_jaffal_results.csv"
        File? output_confident_fusion_fa = "~{sample_id}_jaffal_results.fasta"
        File? output_barcodes_txt = "~{sample_id}_jaffal_results.txt"
        File? output_transcriptome_paf = "~{sample_id}.paf"
        File? output_genome_paf = "~{sample_id}_genome.paf"
        File? output_genome_psl = "~{sample_id}_genome.psl"
        File? output_fusion_fasta = "~{sample_id}_fusions.fa"
        File? output_reads = "~{sample_id}.reads"
        File? output_summary = "~{sample_id}.summary"
        File? output_log = "~{sample_id}_commandlog.txt"
        File? monitoring_log = "~{sample_id}_monitoring.log"
    }
}

task seqtkTask {
    input {
        File? input_fastq
        String sample_id
        String platform

        RuntimeEnvironment runtime_environment
        File monitoring_script 
    }

    command <<<
        bash ~{monitoring_script} > ~{sample_id}_monitoring.log &
        if [[ ~{platform} == 'pacbio' ]]; then
            seqtk seq -A ~{input_fastq}  | awk '{if (/^>/){print $1"-molecule"} else {print}}' | gzip > ~{sample_id}.fa.gz
        else
            seqtk seq -A ~{input_fastq}  | awk '{if (/^>/){print $1"-ccs"} else {print}}' | gzip > ~{sample_id}.fa.gz
        fi
    >>>

    runtime {
        disks: "local-disk ~{runtime_environment.disk_space} HDD"
        memory: "~{runtime_environment.mem} GB"
        cpu: runtime_environment.cpu
        preemptible: runtime_environment.preemptible
        bootDiskSizeGb: runtime_environment.boot_disk_size
        docker: runtime_environment.docker_image
    }

    output {
        File? fa = "~{sample_id}.fa.gz"
        File? monitoring_log = "~{sample_id}_monitoring.log"
    }
}

task pbmm2 {
    input {
        File? input_fasta
        String sample_id
        File pbmmi
        String pbmm2_platform = 'pacbio'

        File monitoring_script
        RuntimeEnvironment runtime_environment
    }

    command <<<
        bash ~{monitoring_script} > ~{sample_id}_monitoring.log &
        if [[ ~{pbmm2_platform} == 'pacbio' ]]; then
            pbmm2 align --num-threads ~{runtime_environment.cpu} --preset ISOSEQ --sort ~{pbmmi} ~{input_fasta} ~{sample_id}.bam
        else
            pbmm2 align --num-threads ~{runtime_environment.cpu} --preset CCS --sort ~{pbmmi} ~{input_fasta} ~{sample_id}.bam
        fi
    >>>

    runtime {
        disks: "local-disk ~{runtime_environment.disk_space} HDD"
        memory: "~{runtime_environment.mem} GB"
        cpu: runtime_environment.cpu
        preemptible: runtime_environment.preemptible
        bootDiskSizeGb: runtime_environment.boot_disk_size
        docker: runtime_environment.docker_image
    }

    output {
        File? bam = "~{sample_id}.bam"
        File? monitoring_log = "~{sample_id}_monitoring.log"
    }
}

task pbfusionTask {
    input {
        File? bam
        File pbfusion_gtf_bin
        String sample_id

        File monitoring_script
        RuntimeEnvironment runtime_environment

        Int threshold=0
        String strand_option=" --strand-insensitive "
    }

    command <<<
        bash ~{monitoring_script} > ~{sample_id}_monitoring.log &

        pbfusion discover --gtf ~{pbfusion_gtf_bin} --output-prefix ~{sample_id}_pbfusion --threads ~{runtime_environment.cpu} -v --min-coverage ~{threshold} ~{bam}

    >>>

    runtime {
        disks: "local-disk ~{runtime_environment.disk_space} HDD"
        memory: "~{runtime_environment.mem} GB"
        cpu: runtime_environment.cpu
        preemptible: runtime_environment.preemptible
        bootDiskSizeGb: runtime_environment.boot_disk_size
        docker: runtime_environment.docker_image
    }

    output {
        File? output_confident_fusion = "~{sample_id}_pbfusion.breakpoints.groups.bed"
        File? monitoring_log = "~{sample_id}_monitoring.log"
    }
}

task CTAT_LR_FUSION_TASK {

  input {
       String sample_name
       File transcripts
       File genome_lib_tar
       Int min_per_id
       Int min_J
       Int min_sumJS    
       Int min_novel_junction_support
       File? illumina_left_fq
       File? illumina_right_fq

       RuntimeEnvironment runtime_environment
       File monitoring_script
  }
  
  command <<<

    set -ex

    bash ~{monitoring_script} > ~{sample_name}_monitoring.log &

    # untar the genome lib
    tar xvf ~{genome_lib_tar}
    # rm ~{genome_lib_tar}
    
    # ctat-LR-fusion

    ctat-LR-fusion --version

    ctat-LR-fusion -T ~{transcripts} \
                --genome_lib_dir ctat_genome_lib_build_dir \
                --min_J ~{min_J}  --min_sumJS ~{min_sumJS} --min_novel_junction_support ~{min_novel_junction_support} \
                --min_per_id ~{min_per_id} \
                --CPU ~{runtime_environment.cpu} \
                --vis \
                ~{"--left_fq " + illumina_left_fq} ~{"--right_fq " + illumina_right_fq } \
                -o ctat_LR_fusion_outdir


    mv ctat_LR_fusion_outdir/ctat-LR-fusion.fusion_predictions.preliminary.tsv ~{sample_name}.ctat-LR-fusion.fusion_predictions.preliminary.tsv
    mv ctat_LR_fusion_outdir/ctat-LR-fusion.fusion_predictions.preliminary.abridged.tsv ~{sample_name}.ctat-LR-fusion.fusion_predictions.preliminary.abridged.tsv
  
    mv ctat_LR_fusion_outdir/ctat-LR-fusion.fusion_predictions.tsv ~{sample_name}.ctat-LR-fusion.fusion_predictions.tsv
    mv ctat_LR_fusion_outdir/ctat-LR-fusion.fusion_predictions.abridged.tsv ~{sample_name}.ctat-LR-fusion.fusion_predictions.abridged.tsv 

    mv ctat_LR_fusion_outdir/ctat-LR-fusion.fusion_inspector_web.html ~{sample_name}.ctat-LR-fusion.fusion_inspector_web.html

    mv ctat_LR_fusion_outdir/fusion_intermediates_dir/IGV_prep/igv.genome.fa ~{sample_name}.ctat-LR-fusion.igv.genome.fa
    mv ctat_LR_fusion_outdir/fusion_intermediates_dir/IGV_prep/igv.genome.fa.fai ~{sample_name}.ctat-LR-fusion.igv.genome.fa.fai
    mv ctat_LR_fusion_outdir/fusion_intermediates_dir/IGV_prep/igv.annot.gtf ~{sample_name}.ctat-LR-fusion.igv.annot.gtf
    mv ctat_LR_fusion_outdir/fusion_intermediates_dir/IGV_prep/igv.annot.bed ~{sample_name}.ctat-LR-fusion.igv.annot.bed
    mv ctat_LR_fusion_outdir/fusion_intermediates_dir/IGV_prep/igv.LR.sorted.bam ~{sample_name}.ctat-LR-fusion.igv.LR.sorted.bam
    mv ctat_LR_fusion_outdir/fusion_intermediates_dir/IGV_prep/igv.LR.sorted.bam.bai ~{sample_name}.ctat-LR-fusion.igv.LR.sorted.bam.bai
    mv ctat_LR_fusion_outdir/fusion_intermediates_dir/IGV_prep/igv.pfam.bed ~{sample_name}.ctat-LR-fusion.igv.pfam.bed
    mv ctat_LR_fusion_outdir/fusion_intermediates_dir/IGV_prep/igv.seqsimilar.bed ~{sample_name}.ctat-LR-fusion.igv.seqsimilar.bed
    mv ctat_LR_fusion_outdir/fusion_intermediates_dir/IGV_prep/igv.LR.breakoint.roi.bed ~{sample_name}.ctat-LR-fusion.igv.LR.breakoint.roi.bed

    tar -zcvhf ~{sample_name}.ctat-LR-fusion.igv.tar.gz ~{sample_name}.ctat-LR-fusion.igv.*
    
    >>>
    
    output {
      File? fusion_report="~{sample_name}.ctat-LR-fusion.fusion_predictions.tsv"
      File? fusion_report_abridged="~{sample_name}.ctat-LR-fusion.fusion_predictions.abridged.tsv"

      File? prelim_fusion_report="~{sample_name}.ctat-LR-fusion.fusion_predictions.preliminary.tsv"
      File? prelim_fusion_report_abridged="~{sample_name}.ctat-LR-fusion.fusion_predictions.preliminary.abridged.tsv"

      File? fusion_report_html="~{sample_name}.ctat-LR-fusion.fusion_inspector_web.html"
      File? igv_tar="~{sample_name}.ctat-LR-fusion.igv.tar.gz"
      File? monitoring_log = "~{sample_name}_monitoring.log"
    }
    

    runtime {
      docker: "~{runtime_environment.docker_image}"
      disks: "local-disk " + runtime_environment.disk_space + " HDD"
      memory: "~{runtime_environment.mem} GB"
      cpu: runtime_environment.cpu
      preemptible: runtime_environment.preemptible
    }
}


task raise_exception {
    input {
        String msg
        RuntimeEnvironment runtime_environment
    }
    command {
        echo -e "\n* Error: ${msg}\n" >&2
        exit 2
    }
    output {
        String error_msg = '${msg}'
    }
    runtime {
        maxRetries : 0
        cpu : 1
        memory : '2 GB'
        time : 4
        disks : 'local-disk 10 SSD'
        docker : runtime_environment.docker_image
    }
}
