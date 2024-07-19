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

    tar xvf ~{arriba_genome_ref_tar}

    /arriba*/run_arriba.sh ArribaReferences/STAR_index_GRCh38_GENCODE38 ArribaReferences/*.gtf ArribaReferences/*.fa  ArribaReferences/blacklist_*.tsv.gz ArribaReferences/known_fusions_*.tsv.gz ArribaReferences/protein_domains_*.gff3 ~{threads} ~{left_fq} ~{right_fq}  

    mv fusions.tsv ~{sample_id}.Arriba.fusions.tsv
      

  >>>

  output {
     File arriba_fusions = "~{sample_id}.Arriba.fusions.tsv"
  }

  

  runtime {
    #preemptible: preemptible
    disks: "local-disk " +  disk_space + " HDD"
    docker: docker
    cpu: threads
    memory: "50G"
  }


}



