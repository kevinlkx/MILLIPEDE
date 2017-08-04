#!/bin/bash

genome_dir="/home/home4/kevinluo/MNaseData/Model/PWM/genomeseq/sacCer2_2008"

## for each chr, extract the genome sequence and save in the database
for chr_num in {1..16}; do
  echo perl extractGenomeSeq.pl \
  --seqfile=${genome_dir}/chr${chr_num}.fa \
  --chrNum=${chr_num}
done

