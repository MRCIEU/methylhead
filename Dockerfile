# Use the official Ubuntu 22.04 image as a base
FROM ubuntu:22.04
ARG DEBIAN_FRONTEND=noninteractive


LABEL Author="Onur Oztornaci"
LABEL Version="1.0"
LABEL Description="This container includes FastQC, Trim Galore, Bismark, Bowtie, Bowtie2, BWA, BWAmeth, Samtools, Sambamba, Picard, Multiqc and MethylDackel for DNAm analysis."


# Install necessary system dependencies
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    unzip \
    bzip2 \
    build-essential \
    zlib1g-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libghc-bzlib-dev \
    libhdf5-dev \
    python3 \
    python3-pip \
    openjdk-8-jdk \
    git \
    libxml2-dev \
    libcairo2-dev \
    libopenblas-dev \
    libfreetype6-dev \
    libglib2.0-dev \
    libpng-dev \
    libjpeg-dev \
    software-properties-common \
    dirmngr \
    gpg-agent \
    sudo \
    gnupg \
    ca-certificates \
    r-base \
    r-base-dev \
    libgsl0-dev \
    libglu1-mesa-dev \
    && rm -rf /var/lib/apt/lists/*


# Install R packages
RUN Rscript -e 'install.packages("BiocManager", repos="https://cran.rstudio.com/");    BiocManager::install(c("methylKit","IlluminaHumanMethylation450kanno.ilmn12.hg19","IlluminaHumanMethylation450kmanifest"),dependencies=TRUE, ask=FALSE)' && \
    Rscript -e "install.packages(c('ggplot2', 'remotes', 'dplyr', 'tidyr', 'data.table'), repos='https://cran.rstudio.com/')" && \
    Rscript -e "remotes::install_github('perishky/meffonym')"

# Install FastQC
RUN wget -q http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip && \
    unzip fastqc_v0.11.9.zip && \
    chmod +x FastQC/fastqc && \
    mv FastQC /opt/fastqc && \
    ln -s /opt/fastqc/fastqc /usr/local/bin/fastqc

# Install Trim Galore
RUN wget -q https://github.com/FelixKrueger/TrimGalore/archive/0.6.6.zip && \
    unzip 0.6.6.zip && \
    mv TrimGalore-0.6.6 /opt/trim_galore && \
    ln -s /opt/trim_galore/trim_galore /usr/local/bin/trim_galore && \
    chmod +x /usr/local/bin/trim_galore

# Install Bismark
RUN wget -q https://github.com/FelixKrueger/Bismark/archive/0.22.3.zip && \
    unzip 0.22.3.zip && \
    mv Bismark-0.22.3 /opt/bismark && \
    ln -s /opt/bismark/bismark /usr/local/bin/bismark && \
    chmod +x /usr/local/bin/bismark

# Install Bowtie
RUN wget -q https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.3.0/bowtie-1.3.0-linux-x86_64.zip && \
    unzip bowtie-1.3.0-linux-x86_64.zip && \
    mv bowtie-1.3.0-linux-x86_64 /opt/bowtie && \
    ln -s /opt/bowtie/bowtie /usr/local/bin/bowtie && \
    chmod +x /usr/local/bin/bowtie

# Install Bowtie2
RUN wget -q https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.2/bowtie2-2.4.2-linux-x86_64.zip && \
    unzip bowtie2-2.4.2-linux-x86_64.zip && \
    mv bowtie2-2.4.2-linux-x86_64 /opt/bowtie2 && \
    ln -s /opt/bowtie2/bowtie2 /usr/local/bin/bowtie2 && \
    chmod +x /usr/local/bin/bowtie2

# Install BWA
RUN git clone https://github.com/lh3/bwa.git && \
    cd bwa && \
    make && \
    mv bwa /opt/bwa && \
    ln -s /opt/bwa /usr/local/bin/bwa && \
    chmod +x /usr/local/bin/bwa && \
    cd ..

# Install toolshed dependency
RUN wget -q https://files.pythonhosted.org/packages/08/bd/d6259c8d882f6783fdfe69bfaed628afb1ddd291f8a1ac6176d27c62860c/toolshed-0.4.6.tar.gz && \
    tar -xzf toolshed-0.4.6.tar.gz && \
    cd toolshed-0.4.6 && \
    python3 setup.py install && \
    cd ..

# Install BWAmeth
RUN wget -q https://github.com/brentp/bwa-meth/archive/master.zip && \
    unzip master.zip && \
    cd bwa-meth-master && \
    python3 setup.py install && \
    cp bwameth.py /usr/local/bin/bwameth && \
    chmod +x /usr/local/bin/bwameth && \
    cd ..

# Install Samtools
RUN wget -q https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2 && \
    tar -xjf samtools-1.10.tar.bz2 && \
    cd samtools-1.10 && \
    ./configure && \
    make && \
    make install && \
    ln -sf /usr/local/bin/samtools /usr/bin/samtools && \
    cd ..

# Install Sambamba
RUN wget -q https://github.com/biod/sambamba/releases/download/v1.0.1/sambamba-1.0.1-linux-amd64-static.gz && \
    gunzip sambamba-1.0.1-linux-amd64-static.gz && \
    chmod +x sambamba-1.0.1-linux-amd64-static && \
    mv sambamba-1.0.1-linux-amd64-static /usr/local/bin/sambamba

# Install htslib dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    libssl-dev \
    tar \
    autoconf \
    automake \
    libtool \
    && rm -rf /var/lib/apt/lists/*

# Install htslib
RUN wget https://github.com/samtools/htslib/releases/download/1.11/htslib-1.11.tar.bz2 && \
    tar -xjf htslib-1.11.tar.bz2 && \
    cd htslib-1.11 && \
    ./configure && \
    make && \
    make install && \
    cd ..

# Install libBigWig
RUN git clone https://github.com/dpryan79/libBigWig.git && \
    cd libBigWig && \
    make && \
    make install && \
    cd ..

# Install MethylDackel from source
RUN git clone https://github.com/dpryan79/MethylDackel.git && \
    cd MethylDackel && \
    make LIBBIGWIG="/usr/local/lib/libBigWig.a" && \
    cp MethylDackel /usr/local/bin/MethylDackel && \
    chmod +x /usr/local/bin/MethylDackel && \
    cd ..

# Set LD_LIBRARY_PATH to include htslib and libBigWig libraries
ENV LD_LIBRARY_PATH="/usr/local/lib:${LD_LIBRARY_PATH}"

# Install Picard
RUN wget -O /opt/picard.jar https://github.com/broadinstitute/picard/releases/download/2.23.8/picard.jar

# Create a script to run Picard directly
RUN echo '#!/bin/bash' > /usr/local/bin/picard && \
    echo 'java -jar /opt/picard.jar "$@"' >> /usr/local/bin/picard && \
    chmod +x /usr/local/bin/picard

# Install Cutadapt
RUN pip3 install cutadapt

# Install MultiQC
RUN pip3 install multiqc

# Set PATH for tools

ENV PATH="/usr/local/bin:/opt/fastqc:/opt/trim_galore:/opt/bismark:/opt/bowtie:/opt/bowtie2:/opt:${PATH}"
RUN mkdir -p /opt/nextflow && chmod 777 /opt/nextflow

CMD ["bash"]

