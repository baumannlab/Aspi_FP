#!/bin/bash
java -Xmx4g -jar /home/dut/bin/picard-tools-1.119//MarkDuplicates.jar I=/home/dut/projects/tigris/heterozygosity/dwn_sample_atig_122/preprocessing/merged_bams/Atig_122.merged.bam O=/home/dut/projects/tigris/heterozygosity/dwn_sample_atig_122/preprocessing/dedup_merged/Atig_122.merged.dedup.bam M=/home/dut/projects/tigris/heterozygosity/dwn_sample_atig_122/preprocessing/dedup_merged/Atig_122.merged_metrics.txt &
proc1=$!
wait "$proc1"
