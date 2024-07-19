version 1.0

workflow Arriba_wf {

  input {
    File arriba_genome_ref_tar = "gs://mdl-refs/GRCh38/ArribaGRCh38/ArribaReferences.tar"
    String docker = "uhrigs/arriba:2.4.0"

    String sample_id
    File left_fq
    File? right_fq
  }

   
  call Arriba_task {
    input:
      arriba_genome_ref_tar = arriba_genome_ref_tar,
      docker = docker,
      sample_id = sample_id,
      left_fq = left_fq,
      right_fq = right_fq
  }

  output {
    File arriba_fusions = Arriba_task.arriba_fusions
  }
}
 

task Arriba_task {
  input {
    File arriba_genome_ref_tar
    String docker

    String sample_id
    File left_fq
    File? right_fq

    Int threads = 8

  }

  Int disk_space = ceil(1+size(arriba_genome_ref_tar, "GB")*10 + size(left_fq, "GB")*2)

  command <<<

    set -ex

    # tar xvf ~{arriba_genome_ref_tar}


    echo /arriba*/run_arriba.sh ArribaReferences/STAR_index_GRCh38_GENCODE38 ArribaReferences/*.gtf ArribaReferences/*.fa  ArribaReferences/blacklist_*.tsv.gz ArribaReferences/known_fusions_*.tsv.gz ArribaReferences/protein_domains_*.gff3 ~{threads} ~{left_fq} ~{right_fq}  



    
    
    #/arriba*/run_arriba.sh  STAR_index_hs37d5viral_GENCODE19/ GENCODE19.gtf hs37d5viral.fa database/blacklist_hg19_hs37d5_GRCh37_v2.4.0.tsv.gz database/known_fusions_hg19_hs37d5_GRCh37_v2.4.0.tsv.gz database/protein_domains_hg19_hs37d5_GRCh37_v2.4.0.gff3 8 test/read1.fastq.gz test/read2.fastq.gz





    # /arriba*/run_arriba.sh /references/STAR_index_* /references/*.gtf /references/*.fa /references/blacklist_*.tsv.gz /references/known_fusions_*.tsv.gz /references/protein_domains_*.gff3 ${THREADS-8} /read1.fastq.gz $(ls /read2.fastq.gz 2> /dev/null)' > /usr/local/bin/arriba.sh && \


  >>>

  output {
     File arriba_fusions = "~{sample_id}.fusions.tsv"
  }
    

  runtime {
    #preemptible: preemptible
    disks: "local-disk " +  disk_space + " HDD"
    docker: docker
    cpu: threads
    memory: "50G"
  }


}



