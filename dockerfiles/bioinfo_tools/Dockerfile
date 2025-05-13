FROM debian:latest

# Update and install system packages
RUN apt update && apt upgrade -y && apt install -y \
    curl unzip pigz pbzip2 default-jre perl build-essential vim-common \
    zlib1g-dev xxd libncurses-dev liblzma-dev libbz2-dev libcurl4-openssl-dev \
    libssl-dev libncurses5-dev libsqlite3-dev libreadline-dev libtk8.6 \
    libgdm-dev libdb4o-cil-dev libpcap-dev jq git

# Download and install Atria
RUN curl -L -O https://github.com/cihga39871/Atria/releases/download/v4.0.3/atria-4.0.3-linux-ubuntu22.tar.gz && \
    tar -zxf atria-4.0.3-linux-ubuntu22.tar.gz && \
    rm atria-4.0.3-linux-ubuntu22.tar.gz
ENV PATH="/atria-4.0.3/bin:$PATH"

# Download and install FastQC
RUN curl -L -O https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip && \
    unzip fastqc_v0.12.1.zip && \
    rm fastqc_v0.12.1.zip

# Download and install featureCounts
RUN curl -L -O https://sourceforge.net/projects/subread/files/subread-2.0.6/subread-2.0.6-Linux-x86_64.tar.gz && \
    tar -zxf subread-2.0.6-Linux-x86_64.tar.gz && \
    rm subread-2.0.6-Linux-x86_64.tar.gz
ENV PATH="/subread-2.0.6-Linux-x86_64/bin:$PATH"

# Download and compile STAR
RUN curl -L -O https://github.com/alexdobin/STAR/archive/2.7.11b.tar.gz && \
    tar -xzf 2.7.11b.tar.gz && \
    rm 2.7.11b.tar.gz && \
    make -C /STAR-2.7.11b/source STAR
ENV PATH="/STAR-2.7.11b/bin/Linux_x86_64_static:$PATH"

# Download and compile samtools
RUN curl -L -O https://github.com/samtools/samtools/releases/download/1.20/samtools-1.20.tar.bz2 && \
    tar -xvjf samtools-1.20.tar.bz2 && \
    rm samtools-1.20.tar.bz2 && \
    make -C /samtools-1.20 && \
    make -C /samtools-1.20 install

# Install Python 3.11
RUN curl -L -O https://www.python.org/ftp/python/3.11.9/Python-3.11.9.tgz && \
    tar xzf Python-3.11.9.tgz && \
    rm Python-3.11.9.tgz && \
    cd Python-3.11.9 && \
    ./configure --enable-optimizations && \
    make -j 4 && \
    make install

# Install bedtools
RUN curl -L -O https://github.com/arq5x/bedtools2/archive/refs/tags/v2.31.1.zip && \
    unzip v2.31.1.zip && \
    rm v2.31.1.zip && \
    make -C /bedtools2-2.31.1
ENV PATH="/bedtools2-2.31.1/bin:$PATH"

# Install Salmon
RUN curl -L -O https://github.com/COMBINE-lab/salmon/releases/download/v1.10.0/salmon-1.10.0_linux_x86_64.tar.gz && \
    tar xzvf salmon-1.10.0_linux_x86_64.tar.gz && \
    rm salmon-1.10.0_linux_x86_64.tar.gz
ENV PATH="/salmon-latest_linux_x86_64/bin:$PATH"

# Install Python packages from requirement file
COPY . .
RUN pip3 install -r requirements.txt
