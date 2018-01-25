#!/bin/bash
(samtools merge -@ 8 - /home/dut/projects/tigris/heterozygosity/dwn_sample_atig_122/preprocessing/dedup_individuals/Atig_122_L21676_HJ2YHBCXX_2.dedup.bam /home/dut/projects/tigris/heterozygosity/dwn_sample_atig_122/preprocessing/dedup_individuals/Atig_122_L21676_HJ2YHBCXX_1.dedup.bam | samtools sort - -m 10G -@ 8 -T /scratch/dut/Atig_122.temp -o /home/dut/projects/tigris/heterozygosity/dwn_sample_atig_122/preprocessing/merged_bams/Atig_122.merged.bam; samtools index /home/dut/projects/tigris/heterozygosity/dwn_sample_atig_122/preprocessing/merged_bams/Atig_122.merged.bam) &
proc1=$!
wait "$proc1"
