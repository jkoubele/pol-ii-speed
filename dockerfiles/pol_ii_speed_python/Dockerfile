FROM python:3.11-bookworm

# Install external tools and system dependencies
RUN apt-get update && apt-get install -y \
    samtools \
    bedtools \
    gzip \
    pigz \
    make \
    gcc \
    zlib1g-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt .
RUN pip3 install --no-cache-dir -r requirements.txt
