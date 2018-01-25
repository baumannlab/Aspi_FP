#!/bin/bash
seqtk sample -s633 /n/analysis/Baumann/rrs/MOLNG-1575/HJ2YHBCXX/s_2_1_CTCAGA.fastq.gz 0.33 > /home/dut/projects/tigris/heterozygosity/dwn_sample_atig_122/preprocessing/down_sample/s_2_1_CTCAGA.0.33ds.fastq; seqtk sample -s633 /n/analysis/Baumann/rrs/MOLNG-1575/HJ2YHBCXX/s_2_2_CTCAGA.fastq.gz 0.33 > /home/dut/projects/tigris/heterozygosity/dwn_sample_atig_122/preprocessing/down_sample/s_2_2_CTCAGA.0.33ds.fastq &
proc1=$!
seqtk sample -s633 /n/analysis/Baumann/rrs/MOLNG-1575/HJ2YHBCXX/s_1_1_CTCAGA.fastq.gz 0.33 > /home/dut/projects/tigris/heterozygosity/dwn_sample_atig_122/preprocessing/down_sample/s_1_1_CTCAGA.0.33ds.fastq; seqtk sample -s633 /n/analysis/Baumann/rrs/MOLNG-1575/HJ2YHBCXX/s_1_2_CTCAGA.fastq.gz 0.33 > /home/dut/projects/tigris/heterozygosity/dwn_sample_atig_122/preprocessing/down_sample/s_1_2_CTCAGA.0.33ds.fastq &
proc2=$!
wait "$proc1" "$proc2"
