#!/bin/bash
samtools index -b /home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/dedup_merged/A_tigris8450.merged.dedup.bam &
proc1=$!
samtools index -b /home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/dedup_merged/Atig_122.merged.dedup.bam &
proc2=$!
samtools index -b /home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/dedup_merged/Atig003.merged.dedup.bam &
proc3=$!
wait "$proc1" "$proc2" "$proc3"
samtools index -b /home/dut/projects/tigris/heterozygosity/gatk6/preprocessing/dedup_merged/Atig001.merged.dedup.bam &
proc4=$!
wait "$proc1" "$proc2" "$proc3" "$proc4"
