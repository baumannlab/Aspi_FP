#!/bin/bash
samtools index -b /home/dut/projects/tigris/heterozygosity/dwn_sample_atig_122/preprocessing/dedup_merged/Atig_122.merged.dedup.bam &
proc1=$!
wait "$proc1"
