#!/bin/bash
echo "gatk_0_prep_ref_2016-04-21"
time ./gatk_0_prep_ref_2016-04-21.sh &> gatk_0_prep_ref_2016-04-21.out
echo "gatk_1_map_to_ref_2016-04-21"
time ./gatk_1_map_to_ref_2016-04-21.sh &> gatk_1_map_to_ref_2016-04-21.out
echo "gatk_2_dedup_individual_2016-04-21"
time ./gatk_2_dedup_individual_2016-04-21.sh &> gatk_2_dedup_individual_2016-04-21.out
echo "gatk_3_merge_bams_2016-04-21"
time ./gatk_3_merge_bams_2016-04-21.sh &> gatk_3_merge_bams_2016-04-21.out
echo "gatk_4_dedup_merged_2016-04-21"
time ./gatk_4_dedup_merged_2016-04-21.sh &> gatk_4_dedup_merged_2016-04-21.out
echo "gatk_5_index_2016-04-21"
time ./gatk_5_index_2016-04-21.sh &> gatk_5_index_2016-04-21.out
echo "gatk_6_realign_bams_2016-04-21"
time ./gatk_6_realign_bams_2016-04-21.sh &> gatk_6_realign_bams_2016-04-21.out
echo "gatk_7_training_2016-04-21"
time ./gatk_7_training_2016-04-21.sh &> gatk_7_training_2016-04-21.out
echo "gatk_8_recalibrate_bams_2016-04-21"
time ./gatk_8_recalibrate_bams_2016-04-21.sh &> gatk_8_recalibrate_bams_2016-04-21.out
echo "gatk_9_final_variant_calling_2016-04-21"
time ./gatk_9_final_variant_calling_2016-04-21.sh &> gatk_9_final_variant_calling_2016-04-21.out
