#!/bin/bash

if ! command -v minimap2 &> /dev/null; then
    echo "Error: minimap2 not installed"
    exit 4
fi

if ! command -v samtools &> /dev/null; then
    wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2
    cd samtools-1.16.1
    ./configure
    make
    make install
    cd ..
fi

if ! command -v bcftools &> /dev/null; then
    wget https://github.com/samtools/bcftools/releases/download/1.16/bcftools-1.16.tar.bz2
    cd bcftools-1.16.1
    ./configure
    make
    make install
    cd ..
fi

if ! command -v tabix &> /dev/null; then
    wget https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2
    cd htslib-1.16
    ./configure
    make
    make install
    cd ..
fi

if ! command -v medaka &> /dev/null; then
    conda install medaka -c conda-forge -c bioconda
fi

if ! command -v longshot &> /dev/null; then
    conda install longshot -c bioconda
fi

pip install pysam matplotlib biopython
