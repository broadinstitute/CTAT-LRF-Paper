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
docker pull trinityctat/ctat_lr_fusion:0.11.1-dev
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

## pbfusion

``` {bash}
docker pull trinityctat/pbfusion:v0.3.1
```

## FusionSeeker

``` {bash}
# FusionSeeker (v1.0.1, commit 5710dc4)
docker pull qianqin/fusion_seeker
```

# WDL script

We need to upload the `*wdl` on a Terra workspace to process the dataset.
