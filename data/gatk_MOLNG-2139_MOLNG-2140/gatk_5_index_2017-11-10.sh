#!/bin/bash
samtools index -b /n/projects/dut/a_marmorata/parthenogen_heterozygosity/data/gatk_MOLNG-2139_MOLNG-2140/preprocessing/dedup_merged/Atig_4278.merged.dedup.bam &
proc1=$!
samtools index -b /n/projects/dut/a_marmorata/parthenogen_heterozygosity/data/gatk_MOLNG-2139_MOLNG-2140/preprocessing/dedup_merged/A.tig_12512.merged.dedup.bam &
proc2=$!
samtools index -b /n/projects/dut/a_marmorata/parthenogen_heterozygosity/data/gatk_MOLNG-2139_MOLNG-2140/preprocessing/dedup_merged/A.tig_12513.merged.dedup.bam &
proc3=$!
wait "$proc1" "$proc2" "$proc3"
samtools index -b /n/projects/dut/a_marmorata/parthenogen_heterozygosity/data/gatk_MOLNG-2139_MOLNG-2140/preprocessing/dedup_merged/A.tig_9721.merged.dedup.bam &
proc4=$!
samtools index -b /n/projects/dut/a_marmorata/parthenogen_heterozygosity/data/gatk_MOLNG-2139_MOLNG-2140/preprocessing/dedup_merged/Atig_6993.merged.dedup.bam &
proc5=$!
samtools index -b /n/projects/dut/a_marmorata/parthenogen_heterozygosity/data/gatk_MOLNG-2139_MOLNG-2140/preprocessing/dedup_merged/Atig_9177.merged.dedup.bam &
proc6=$!
wait "$proc4" "$proc5" "$proc6"
wait "$proc1" "$proc2" "$proc3" "$proc4" "$proc5" "$proc6"
