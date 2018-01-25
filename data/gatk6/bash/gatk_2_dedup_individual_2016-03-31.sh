#!/bin/bash
java -Xmx4g -jar /home/dut/bin/picard-tools-1.119//MarkDuplicates.jar I=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/map_to_ref_output/A_tigris8450_L13136_HF7YWADXXb_2.bam O=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/dedup_individuals/A_tigris8450_L13136_HF7YWADXXb_2.dedup.bam M=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/dedup_individuals/A_tigris8450_L13136_HF7YWADXXb_2_metrics.txt &
proc1=$!
java -Xmx4g -jar /home/dut/bin/picard-tools-1.119//MarkDuplicates.jar I=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/map_to_ref_output/A_tigris8450_L13136-1_HG7HJADXXb_2.bam O=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/dedup_individuals/A_tigris8450_L13136-1_HG7HJADXXb_2.dedup.bam M=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/dedup_individuals/A_tigris8450_L13136-1_HG7HJADXXb_2_metrics.txt &
proc2=$!
java -Xmx4g -jar /home/dut/bin/picard-tools-1.119//MarkDuplicates.jar I=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/map_to_ref_output/A_tigris8450_L13136-1_HG7HJADXXb_1.bam O=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/dedup_individuals/A_tigris8450_L13136-1_HG7HJADXXb_1.dedup.bam M=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/dedup_individuals/A_tigris8450_L13136-1_HG7HJADXXb_1_metrics.txt &
proc3=$!
wait "$proc1" "$proc2" "$proc3"
java -Xmx4g -jar /home/dut/bin/picard-tools-1.119//MarkDuplicates.jar I=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/map_to_ref_output/A_tigris8450_L13136_HBCBWADXXb_2.bam O=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/dedup_individuals/A_tigris8450_L13136_HBCBWADXXb_2.dedup.bam M=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/dedup_individuals/A_tigris8450_L13136_HBCBWADXXb_2_metrics.txt &
proc4=$!
java -Xmx4g -jar /home/dut/bin/picard-tools-1.119//MarkDuplicates.jar I=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/map_to_ref_output/A_tigris8450_L13136_HBCBWADXXb_1.bam O=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/dedup_individuals/A_tigris8450_L13136_HBCBWADXXb_1.dedup.bam M=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/dedup_individuals/A_tigris8450_L13136_HBCBWADXXb_1_metrics.txt &
proc5=$!
java -Xmx4g -jar /home/dut/bin/picard-tools-1.119//MarkDuplicates.jar I=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/map_to_ref_output/A_tigris8450_L13136_HF7YWADXXb_1.bam O=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/dedup_individuals/A_tigris8450_L13136_HF7YWADXXb_1.dedup.bam M=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/dedup_individuals/A_tigris8450_L13136_HF7YWADXXb_1_metrics.txt &
proc6=$!
wait "$proc4" "$proc5" "$proc6"
java -Xmx4g -jar /home/dut/bin/picard-tools-1.119//MarkDuplicates.jar I=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/map_to_ref_output/Atig_122_L21676_HJ2YHBCXX_2.bam O=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/dedup_individuals/Atig_122_L21676_HJ2YHBCXX_2.dedup.bam M=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/dedup_individuals/Atig_122_L21676_HJ2YHBCXX_2_metrics.txt &
proc7=$!
java -Xmx4g -jar /home/dut/bin/picard-tools-1.119//MarkDuplicates.jar I=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/map_to_ref_output/Atig_122_L21676_HJ2YHBCXX_1.bam O=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/dedup_individuals/Atig_122_L21676_HJ2YHBCXX_1.dedup.bam M=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/dedup_individuals/Atig_122_L21676_HJ2YHBCXX_1_metrics.txt &
proc8=$!
java -Xmx4g -jar /home/dut/bin/picard-tools-1.119//MarkDuplicates.jar I=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/map_to_ref_output/Atig003_L13088_HF7YWADXXb_1.bam O=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/dedup_individuals/Atig003_L13088_HF7YWADXXb_1.dedup.bam M=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/dedup_individuals/Atig003_L13088_HF7YWADXXb_1_metrics.txt &
proc9=$!
wait "$proc7" "$proc8" "$proc9"
java -Xmx4g -jar /home/dut/bin/picard-tools-1.119//MarkDuplicates.jar I=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/map_to_ref_output/Atig003_L13088_HBCBWADXXb_1.bam O=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/dedup_individuals/Atig003_L13088_HBCBWADXXb_1.dedup.bam M=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/dedup_individuals/Atig003_L13088_HBCBWADXXb_1_metrics.txt &
proc10=$!
java -Xmx4g -jar /home/dut/bin/picard-tools-1.119//MarkDuplicates.jar I=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/map_to_ref_output/Atig003_L13088_HBCBWADXXb_2.bam O=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/dedup_individuals/Atig003_L13088_HBCBWADXXb_2.dedup.bam M=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/dedup_individuals/Atig003_L13088_HBCBWADXXb_2_metrics.txt &
proc11=$!
java -Xmx4g -jar /home/dut/bin/picard-tools-1.119//MarkDuplicates.jar I=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/map_to_ref_output/Atig003_L13088-1_HG7HJADXXb_1.bam O=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/dedup_individuals/Atig003_L13088-1_HG7HJADXXb_1.dedup.bam M=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/dedup_individuals/Atig003_L13088-1_HG7HJADXXb_1_metrics.txt &
proc12=$!
wait "$proc10" "$proc11" "$proc12"
java -Xmx4g -jar /home/dut/bin/picard-tools-1.119//MarkDuplicates.jar I=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/map_to_ref_output/Atig003_L13088-1_HG7HJADXXb_2.bam O=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/dedup_individuals/Atig003_L13088-1_HG7HJADXXb_2.dedup.bam M=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/dedup_individuals/Atig003_L13088-1_HG7HJADXXb_2_metrics.txt &
proc13=$!
java -Xmx4g -jar /home/dut/bin/picard-tools-1.119//MarkDuplicates.jar I=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/map_to_ref_output/Atig003_L13088_HF7YWADXXb_2.bam O=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/dedup_individuals/Atig003_L13088_HF7YWADXXb_2.dedup.bam M=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/dedup_individuals/Atig003_L13088_HF7YWADXXb_2_metrics.txt &
proc14=$!
java -Xmx4g -jar /home/dut/bin/picard-tools-1.119//MarkDuplicates.jar I=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/map_to_ref_output/Atig001_L13087_HG7HJADXXb_2.bam O=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/dedup_individuals/Atig001_L13087_HG7HJADXXb_2.dedup.bam M=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/dedup_individuals/Atig001_L13087_HG7HJADXXb_2_metrics.txt &
proc15=$!
wait "$proc13" "$proc14" "$proc15"
java -Xmx4g -jar /home/dut/bin/picard-tools-1.119//MarkDuplicates.jar I=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/map_to_ref_output/Atig001_L13087_HG7HJADXXb_1.bam O=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/dedup_individuals/Atig001_L13087_HG7HJADXXb_1.dedup.bam M=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/dedup_individuals/Atig001_L13087_HG7HJADXXb_1_metrics.txt &
proc16=$!
java -Xmx4g -jar /home/dut/bin/picard-tools-1.119//MarkDuplicates.jar I=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/map_to_ref_output/Atig001_L13087_HF7YWADXXb_2.bam O=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/dedup_individuals/Atig001_L13087_HF7YWADXXb_2.dedup.bam M=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/dedup_individuals/Atig001_L13087_HF7YWADXXb_2_metrics.txt &
proc17=$!
java -Xmx4g -jar /home/dut/bin/picard-tools-1.119//MarkDuplicates.jar I=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/map_to_ref_output/Atig001_L13087_HBCBWADXXb_2.bam O=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/dedup_individuals/Atig001_L13087_HBCBWADXXb_2.dedup.bam M=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/dedup_individuals/Atig001_L13087_HBCBWADXXb_2_metrics.txt &
proc18=$!
wait "$proc16" "$proc17" "$proc18"
java -Xmx4g -jar /home/dut/bin/picard-tools-1.119//MarkDuplicates.jar I=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/map_to_ref_output/Atig001_L13087_HBCBWADXXb_1.bam O=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/dedup_individuals/Atig001_L13087_HBCBWADXXb_1.dedup.bam M=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/dedup_individuals/Atig001_L13087_HBCBWADXXb_1_metrics.txt &
proc19=$!
java -Xmx4g -jar /home/dut/bin/picard-tools-1.119//MarkDuplicates.jar I=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/map_to_ref_output/Atig001_L13087_HF7YWADXXb_1.bam O=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/dedup_individuals/Atig001_L13087_HF7YWADXXb_1.dedup.bam M=/home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/dedup_individuals/Atig001_L13087_HF7YWADXXb_1_metrics.txt &
proc20=$!
