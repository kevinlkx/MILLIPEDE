#!/bin/bash

echo There are $# arguments to $0: $*
echo motif_id: $1
echo tfName: $2
echo PWM score: $3

## scan motif matches and save in mySQL database
perl /home/home4/kevinluo/MNaseData/Model/PWM/motif_scan/Perl_Functions/motif_scan_cs.pl \
--n=4 \
--pwm=/home/home4/kevinluo/MNaseData/Model/PWM/motif_scan/motif_files/tf_coord_updated_macisaac/$2.pwm \
--motif_id=$1 \
--desc="4th order $2 PWM with a min score of 2" \
--bfile=/home/home4/kevinluo/MNaseData/Model/PWM/motif_scan/motif_files/bg_all.tab \
--min_score=2

## extract motif matches from the mySQL database
perl /home/home4/kevinluo/MNaseData/Model/PWM/motif_scan/Perl_Functions/get_pwm_dbi_cs.pl \
--motif_id=$1 \
--score=$3

## convert the results into a GenomicRanges object
Rscript /home/home4/kevinluo/MNaseData/Model/PWM/motif_scan/Perl_Functions/pwm2gr.R $1 $2 $3 \
