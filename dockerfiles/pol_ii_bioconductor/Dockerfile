FROM bioconductor/bioconductor_docker:latest

RUN R -e "install.packages(c('tidyverse', 'rjson', 'argparse', 'arrow'))"
RUN R -e  "BiocManager::install('rtracklayer')"

RUN apt-get update && apt-get install -y pigz

# Install bedtools
RUN curl -L -O https://github.com/arq5x/bedtools2/archive/refs/tags/v2.31.1.zip && \
    unzip v2.31.1.zip && \
    rm v2.31.1.zip && \
    make -C /bedtools2-2.31.1
ENV PATH="/bedtools2-2.31.1/bin:$PATH"

RUN R -e  "BiocManager::install('txdbmaker')"
RUN R -e  "BiocManager::install('DESeq2')"