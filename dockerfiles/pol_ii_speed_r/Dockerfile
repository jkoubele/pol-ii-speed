FROM bioconductor/bioconductor_docker:RELEASE_3_21-R-4.5.1

RUN apt-get update && apt-get install -y pigz curl unzip make g++

RUN R -e "install.packages(c('tidyverse', 'rjson', 'argparse', 'arrow', 'ashr'))"
RUN R -e  "BiocManager::install(c('rtracklayer','DESeq2', 'txdbmaker', 'GenomicFeatures', 'tximport'))"

# Install bedtools
RUN curl -L -O https://github.com/arq5x/bedtools2/archive/refs/tags/v2.31.1.zip && \
    unzip v2.31.1.zip && \
    rm v2.31.1.zip && \
    make -C /bedtools2-2.31.1
ENV PATH="/bedtools2-2.31.1/bin:$PATH"
