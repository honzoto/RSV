FROM ubuntu
MAINTAINER Honzo <hz_juanito18@outlook.com>

LABEL description="Basic workflow for RNA-seq Variant Calling"

ARG DEBIAN_FRONTEND=noninteractive
#WORKDIR /app

# Get ubuntu pre-requisites
RUN apt-get update -y
RUN apt-get install -y \
        build-essential \
        bc \
        unzip \
        zip \
        python3 \
        pip \
        python-is-python3 \
        default-jdk \
        default-jre \
        git \
        curl \
        wget \
        curl \
        liblzma-dev \
        libncurses-dev \
        libbz2-dev \
        libcurl3-dev \
        libhdf5-serial-dev \
        libssl-dev \
        libboost-all-dev \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    apt-get clean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/{apt,dpkg,cache,log}/

#ln -s /usr/bin/python3 /usr/bin/python

# Get pre-requisite tools
COPY /opt /opt

# install BBtools
WORKDIR /opt
RUN tar -xvzf BBMap_38.98.tar.gz

# install SDKMan
WORKDIR /opt
RUN curl -s "https://get.sdkman.io" | bash && \
    bash /root/.sdkman/bin/sdkman-init.sh

# install Picard tools
RUN git clone https://github.com/broadinstitute/picard.git
WORKDIR /opt/picard
RUN ./gradlew shadowJar
#RUN java -jar /opt/picard/build/libs/picard.jar

# installing GATK
WORKDIR /opt
RUN unzip gatk-4.2.6.1.zip && \
    cd gatk-4.2.6.1 && \
    bash gatk-completion.sh

# install FastQC
WORKDIR /opt
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip && \
    unzip fastqc_v0.11.9.zip && \
    mv FastQC FastQC-0.11.9 && \
    cd FastQC-0.11.9 && \
    chmod 775 fastqc
ENV PATH /opt/FastQC-0.11.9:$PATH

# install HTSLIB tools (bcftools, samtools)
WORKDIR /opt
RUN tar -xf bcftools-1.16.tar.bz2 && \
    cd bcftools-1.16 && \
    make
ENV PATH /opt/bcftools-1.16/bcftools:$PATH

RUN tar -xf samtools-1.16.1.tar.bz2 && \
    cd samtools-1.16.1 && \
    make
ENV PATH /opt/samtools-1.16.1/samtools:$PATH

# install BWA
WORKDIR /opt
RUN git clone https://github.com/lh3/bwa.git && \
    cd bwa && \
    make
ENV PATH /opt/bwa/bwa:$PATH

# install python packages indicated in requirements.txt
COPY requirements.txt requirements.txt
RUN pip install -r requirements.txt

# Get genome references (we've already indexed them)
COPY /data /data
#WORKDIR /data/INDEX_HG38
#RUN /opt/bwa/bwa index GRCh38.p13_genomic.fna

# CLEANING UP
WORKDIR /opt
RUN rm *.bz2 *.zip *.gz
# is this the right way to add a directory to the path?
ENV PATH /app:$PATH
#RUN alias python='/usr/bin/python3'

# Get fastq files (if applicable)
COPY /fastq /fastq

# Get pipeline and scripts
COPY /app /app
WORKDIR /app





