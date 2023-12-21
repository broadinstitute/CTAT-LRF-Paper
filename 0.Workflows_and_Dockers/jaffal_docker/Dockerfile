FROM ubuntu:18.04
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y \
    && apt-get install -y \
        bowtie2 \
        bzip2 \
        g++ \
        git \
        gzip \
        libncurses5-dev \
        libpng-dev \
        libtool \
        lmod \
        make \
        openjdk-8-jdk \
        python \
        r-base \
        r-base-dev \
        time \
        trimmomatic \
        unzip \
        wget \
        zlib1g-dev \
    && rm -rf /tmp/downloaded_packages/ /tmp/*.rds \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt

# Set Standard settings
ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64
ENV CLASSPATH .
ENV CP ${CLASSPATH}
ENV BASH_ENV /usr/share/lmod/lmod/init/bash
ENV PATH /usr/bin:/bin:/usr/local/bin:/opt/bin:/opt/bbmap
#ENV MODULEPATH

# install R dependencies (required by JAFFA)
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R")' -e 'biocLite("IRanges")'

# install jaffa
RUN git clone https://github.com/Oshlack/JAFFA.git -b master
RUN mkdir -p JAFFA/tools/bin
RUN g++ -std=c++11 -O3 -o JAFFA/tools/bin/process_transcriptome_align_table JAFFA/src/process_transcriptome_align_table.c++
RUN g++ -O3 -o JAFFA/tools/bin/extract_seq_from_fasta JAFFA/src/extract_seq_from_fasta.c++
RUN g++ -std=c++11 -O3 -o JAFFA/tools/bin/make_simple_read_table JAFFA/src/make_simple_read_table.c++
#RUN g++ -std=c++11 -O3 -o JAFFA/bin/bypass_genomic_alignment JAFFA/src/bypass_genomic_alignment.c++
ENV PATH ${PATH}:/opt/JAFFA/tools/bin
# set the tools
RUN cd JAFFA; bash install_linux64.sh
RUN chmod 644 JAFFA/tools.groovy

RUN cd JAFFA; wget -c https://figshare.com/ndownloader/files/25410494 -O JAFFA_REFERENCE_FILES_HG38_GENCODE22.V2.tar.gz; tar -zxvf JAFFA_REFERENCE_FILES_HG38_GENCODE22.V2.tar.gz; rm JAFFA_REFERENCE_FILES_HG38_GENCODE22.V2.tar.gz

