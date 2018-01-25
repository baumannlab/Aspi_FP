#!/bin/bash
ln -s /home/dut/projects/tigris/genome_annotation/fasta/tigris_scaffolds_filt_10000.fa /home/dut/projects/tigris/heterozygosity/dwn_sample_atig_122/reference/
bwa index -a bwtsw /home/dut/projects/tigris/heterozygosity/dwn_sample_atig_122/reference/tigris_scaffolds_filt_10000.fa
samtools faidx /home/dut/projects/tigris/heterozygosity/dwn_sample_atig_122/reference/tigris_scaffolds_filt_10000.fa
java -jar /home/dut/bin/picard-tools-1.119/CreateSequenceDictionary.jar REFERENCE=/home/dut/projects/tigris/heterozygosity/dwn_sample_atig_122/reference/tigris_scaffolds_filt_10000.fa OUTPUT=/home/dut/projects/tigris/heterozygosity/dwn_sample_atig_122/reference/tigris_scaffolds_filt_10000.dict
