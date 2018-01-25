#!/bin/bash
java -Xmx4g -jar /home/dut/bin/picard-tools-1.119//MarkDuplicates.jar I=/home/dut/projects/tigris/heterozygosity/dwn_sample_atig_122/preprocessing/map_to_ref_output/Atig_122_L21676_HJ2YHBCXX_2.bam O=/home/dut/projects/tigris/heterozygosity/dwn_sample_atig_122/preprocessing/dedup_individuals/Atig_122_L21676_HJ2YHBCXX_2.dedup.bam M=/home/dut/projects/tigris/heterozygosity/dwn_sample_atig_122/preprocessing/dedup_individuals/Atig_122_L21676_HJ2YHBCXX_2_metrics.txt &
proc1=$!
java -Xmx4g -jar /home/dut/bin/picard-tools-1.119//MarkDuplicates.jar I=/home/dut/projects/tigris/heterozygosity/dwn_sample_atig_122/preprocessing/map_to_ref_output/Atig_122_L21676_HJ2YHBCXX_1.bam O=/home/dut/projects/tigris/heterozygosity/dwn_sample_atig_122/preprocessing/dedup_individuals/Atig_122_L21676_HJ2YHBCXX_1.dedup.bam M=/home/dut/projects/tigris/heterozygosity/dwn_sample_atig_122/preprocessing/dedup_individuals/Atig_122_L21676_HJ2YHBCXX_1_metrics.txt &
proc2=$!
wait "$proc1" "$proc2"
