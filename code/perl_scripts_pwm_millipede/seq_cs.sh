#!/bin/bash

for i in {1..16}
do
perl /home/home4/kevinluo/MNaseData/Model/PWM/genomeseq/extractGenomeSeq.pl \
--seqfile="/home/home4/kevinluo/MNaseData/Model/PWM/genomeseq/sacCer2_2008/chr"$i".fa" \
--chrNum=$i
done

