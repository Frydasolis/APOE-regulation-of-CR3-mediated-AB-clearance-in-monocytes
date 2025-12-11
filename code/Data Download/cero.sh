#!/bin/bash
mkdir -p /datos/lymphocytes/monocytes/ad/fastqfiles
cd /datos/lymphocytes/monocytes/ad/fastqfiles || exit

echo "Starting download of AD samples..."
while read accession; do
    echo "Downloading $accession..."
    fasterq-dump $accession
done < AD_samples.txt

echo "Starting download of Control samples..."
while read accession; do
    echo "Downloading $accession..."
    fasterq-dump $accession
done < Control_samples.txt

echo "Compressing FASTQ files..."
gzip *.fastq

echo "All downloads and compression complete."

