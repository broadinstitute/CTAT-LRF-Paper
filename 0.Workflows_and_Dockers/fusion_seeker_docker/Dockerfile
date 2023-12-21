FROM ubuntu

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get -qq update

RUN apt-get install -qq -y wget git build-essential zlib1g-dev curl gcc g++ make

ENV LANG C.UTF-8
ENV SHELL /bin/bash
ENV SRC=/usr/local/src
ENV BIN=/usr/local/bin

RUN /bin/bash -c "curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-\$(uname -m).sh > mambaforge.sh && \
    bash mambaforge.sh -b -p /opt/conda && \
    rm mambaforge.sh"

ENV PATH="/opt/conda/bin:/bsalign:/FusionSeeker:${PATH}"
ARG PATH="/opt/conda/bin:/bsalign:/FusionSeeker:${PATH}"

RUN mamba create -q -y -c bioconda -n fusionseeker minimap2 pysam

RUN git clone https://github.com/ruanjue/bsalign.git && cd /bsalign && make
RUN git clone https://github.com/Maggi-Chen/FusionSeeker.git

ENV PATH="/opt/conda/envs/fusionseeker/bin:/bsalign:/FusionSeeker:${PATH}"

RUN apt-get install -qq -y samtools 

RUN echo "source activate fusionseeker" > ~/.bashrc
