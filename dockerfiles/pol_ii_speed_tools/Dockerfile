FROM debian:bookworm

# Update and install system packages
RUN apt update && apt upgrade -y && apt install -y \
    curl unzip pigz pbzip2 default-jre perl build-essential vim-common \
    zlib1g-dev xxd libncurses-dev liblzma-dev libbz2-dev libcurl4-openssl-dev \
    libssl-dev libncurses5-dev libsqlite3-dev libreadline-dev libtk8.6 \
    libgdm-dev libdb4o-cil-dev libpcap-dev jq git cmake libgsl-dev

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

# Install and compile MashMap
RUN curl -L -O https://github.com/marbl/MashMap/archive/refs/tags/v3.1.3.tar.gz && \
    tar xzvf v3.1.3.tar.gz && \
    rm v3.1.3.tar.gz && \
    cmake -S MashMap-3.1.3 -B MashMap-3.1.3/build -DCMAKE_BUILD_TYPE=Release && \
    cmake --build MashMap-3.1.3/build
ENV PATH="/MashMap-3.1.3/build/bin:$PATH"
