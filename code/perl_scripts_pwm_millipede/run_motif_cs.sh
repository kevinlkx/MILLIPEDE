#!/bin/bash

# echo There are $# arguments to $0: $*

motif_id=$1
tf_name=$2
pwm_score=$3

echo motif_id: ${motif_id}
echo tf_name: ${tf_name}
echo PWM_score: ${pwm_score}

## scan motif matches and save in mySQL database
echo perl /home/home4/kevinluo/MNaseData/Model/PWM/motif_scan/Perl_Functions/motif_scan_cs.pl \
--n=4 \
--pwm=/home/home4/kevinluo/MNaseData/Model/PWM/motif_scan/motif_files/tf_coord_updated_macisaac/${tf_name}.pwm \
--motif_id=${motif_id} \
--desc="4th order ${tf_name} PWM with a min score of 2" \
--bfile=/home/home4/kevinluo/MNaseData/Model/PWM/motif_scan/motif_files/bg_all.tab \
--min_score=2

## extract motif matches from the mySQL database
echo perl /home/home4/kevinluo/MNaseData/Model/PWM/motif_scan/Perl_Functions/get_pwm_dbi_cs.pl \
--motif_id=${motif_id} \
--score=${pwm_score}

## convert the results into a GenomicRanges object
echo Rscript /home/home4/kevinluo/MNaseData/Model/PWM/motif_scan/Perl_Functions/pwm2gr.R ${motif_id} ${tf_name} ${pwm_score}


