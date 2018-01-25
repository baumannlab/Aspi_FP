#!/bin/bash
java -Xmx4g -jar /home/dut/bin/picard-tools-1.119//MarkDuplicates.jar I=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/merged_bams/A_tigris8450.merged.bam O=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/dedup_individuals/A_tigris8450.merged.dedup.bam M=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/dedup_individuals/A_tigris8450.merged_metrics.txt &
proc1=$!
java -Xmx4g -jar /home/dut/bin/picard-tools-1.119//MarkDuplicates.jar I=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/merged_bams/Atig_122.merged.bam O=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/dedup_individuals/Atig_122.merged.dedup.bam M=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/dedup_individuals/Atig_122.merged_metrics.txt &
proc2=$!
java -Xmx4g -jar /home/dut/bin/picard-tools-1.119//MarkDuplicates.jar I=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/merged_bams/Atig003.merged.bam O=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/dedup_individuals/Atig003.merged.dedup.bam M=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/dedup_individuals/Atig003.merged_metrics.txt &
proc3=$!
wait "$proc1" "$proc2" "$proc3"
java -Xmx4g -jar /home/dut/bin/picard-tools-1.119//MarkDuplicates.jar I=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/merged_bams/Atig001.merged.bam O=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/dedup_individuals/Atig001.merged.dedup.bam M=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/dedup_individuals/Atig001.merged_metrics.txt &
proc4=$!
wait "$proc1" "$proc2" "$proc3" "$proc4"
