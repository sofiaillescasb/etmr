#!/bin/bash

for srr in $(cat SRR_Acc_List.txt); do
echo "Downloading $srr..."
fasterq-dump $srr --split-files --threads 8
gzip ${srr}_1.fastq ${srr}_2.fastq
echo "Finished $srr"
done
