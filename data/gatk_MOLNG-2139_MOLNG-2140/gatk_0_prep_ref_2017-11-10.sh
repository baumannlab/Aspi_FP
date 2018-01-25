#!/bin/bash
ln -s /home/dut/projects/tigris/genome_annotation/fasta/tigris_scaffolds_filt_10000.fa /n/projects/dut/a_marmorata/parthenogen_heterozygosity/data/gatk_MOLNG-2139_MOLNG-2140/reference/
bwa index -a bwtsw /n/projects/dut/a_marmorata/parthenogen_heterozygosity/data/gatk_MOLNG-2139_MOLNG-2140/reference/tigris_scaffolds_filt_10000.fa
samtools faidx /n/projects/dut/a_marmorata/parthenogen_heterozygosity/data/gatk_MOLNG-2139_MOLNG-2140/reference/tigris_scaffolds_filt_10000.fa
java -jar /home/dut/bin/picard-tools-1.119/CreateSequenceDictionary.jar REFERENCE=/n/projects/dut/a_marmorata/parthenogen_heterozygosity/data/gatk_MOLNG-2139_MOLNG-2140/reference/tigris_scaffolds_filt_10000.fa OUTPUT=/n/projects/dut/a_marmorata/parthenogen_heterozygosity/data/gatk_MOLNG-2139_MOLNG-2140/reference/tigris_scaffolds_filt_10000.dict
