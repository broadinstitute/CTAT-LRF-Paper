Docker images and WDLs for pipeline execution are included for long-read fusion detection software as leveraged here.

# Build from the docker file

Please see individual folder and build with following command:

For example,
``` {bash}
cd jaffal_docker/
docker build -t jaffal .
```

# Pre-Built docker images

## CTAT-LR-fusion

```{bash}
docker pull trinityctat/ctat_lr_fusion:0.13.0
```

## JAFFAL

```{bash}
# version 2.3
docker pull trinityctat/jaffal
```

## LongGF

```{bash}
# version 0.1.2
docker pull trinityctat/longgf
```

LongGF reference file:

```sh
gs://fc-d438f734-542a-445f-96e1-e07b4a7cddbc/references/GRCh38.primary_assembly.genome.fa
gs://fc-d438f734-542a-445f-96e1-e07b4a7cddbc/references/gencode.v22.chr_patch_hapl_scaff.annotation.gtf
```

## pbfusion

``` {bash}
docker pull trinityctat/pbfusion:v0.4.0
```

pbfusion v0.4.0 reference file link:

```sh
gs://fc-d438f734-542a-445f-96e1-e07b4a7cddbc/references/broad-gencode-v22.v0.4.0.gtf.bin
gs://fc-d438f734-542a-445f-96e1-e07b4a7cddbc/references/ref_genome.mmi
```

## FusionSeeker

``` {bash}
# FusionSeeker (v1.0.1, commit 5710dc4)
docker pull qianqin/fusion_seeker
```

FusionSeeker reference file:

```sh
gs://fc-d438f734-542a-445f-96e1-e07b4a7cddbc/references/GRCh38.primary_assembly.genome.fa
gs://fc-d438f734-542a-445f-96e1-e07b4a7cddbc/references/gencode.v22.chr_patch_hapl_scaff.annotation.gtf
```

# WDL script

We need to upload the `*wdl` on a Terra workspace to process the dataset.
