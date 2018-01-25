#!/usr/bin/env python
#Author: Duncan Tormey
#Email: dut@stowers.org or duncantormey@gmail.com

import os
import sys
import datetime
import argparse
from collections import defaultdict
sys.path = ["/home/dut/bin/python_scripts/"] + sys.path
import lims_classes as lims


def make_script_fh(base_name):
    date = datetime.datetime.today()
    date = str(date)
    date = date[:10]
    name = './%s_%s.sh' % (base_name, date)
    fh = open(name, "w")
    return fh, name


def make_output_dir(directory_name):
    working_directory = os.getcwd()
    output_directory = working_directory + "/" + directory_name
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    return output_directory


def make_proc_capture(i):
    proc_capture = 'proc%s=$!\n' % str(i)
    proc = '\"$proc%s\"' % str(i)
    return proc_capture, proc


def prep_reference_commands(fasta_path,  picard_tools_path, write_prep_ref):
    output_dir = make_output_dir("reference")
    fasta_name = os.path.basename(fasta_path)
    fasta_ref_path = output_dir + "/" + fasta_name
    fasta_dict_path = output_dir + "/" + fasta_name.replace(".fa", ".dict")
    create_seq_dict = picard_tools_path + 'CreateSequenceDictionary.jar'
    sym_link = "ln -s %s %s/\n" % (fasta_path, output_dir)
    bwa_index = "bwa index -a bwtsw %s\n" % fasta_ref_path
    samtools_faidx = "samtools faidx %s\n" % fasta_ref_path
    sequence_dir = "java -jar %s REFERENCE=%s OUTPUT=%s\n" % (
        create_seq_dict, fasta_ref_path, fasta_dict_path)

    if write_prep_ref:
        fh, scriptname = make_script_fh("gatk_0_prep_ref")
        fh.write("#!/bin/bash\n")
        fh.write(sym_link)
        fh.write(bwa_index)
        fh.write(samtools_faidx)
        fh.write(sequence_dir)

        fh.close()

    return fasta_ref_path, scriptname


def ret_down_sample_read_pairs(pairs, sample, dwn_frac,outdir):
    #print(pairs)
    dwn_frac = float(dwn_frac)
    r_pair = pairs[0].path
    f_pair = pairs[1].path
    seed=sum([ord(a) for a in sample])
    pairs[0].path = outdir + '/' + os.path.basename(r_pair).replace('.fq.gz', '.%sds.fq'%str(dwn_frac)).replace('.fastq.gz', '.%sds.fastq'%str(dwn_frac))
    pairs[1].path = outdir + '/' + os.path.basename(f_pair).replace('.fq.gz', '.%sds.fq'%str(dwn_frac)).replace('.fastq.gz', '.%sds.fastq'%str(dwn_frac))
    dwn_sample_cmd = 'seqtk sample -s%s %s %s > %s; seqtk sample -s%s %s %s > %s &\n' % (seed, r_pair, dwn_frac, pairs[0].path,
                                                                                         seed, f_pair, dwn_frac, pairs[1].path)
    return dwn_sample_cmd, pairs

def down_sample_cmds(lane_pairs, down_sample, write, write_prefix, outdir):
    print('\n\n')
    print(lane_pairs)
    print('\n\n')
    outdir = make_output_dir(outdir)
    down_samples = dict(down_sample)  
    if write:
        fh, scriptname = make_script_fh(write_prefix)
        fh.write("#!/bin/bash\n")
    n = 1
    all_procs = []
    for sample in lane_pairs:
        if sample in down_samples:
            for  lane, pairs in lane_pairs[sample].items():
                print(sample)
                print(lane)
                print(pairs)
                cmd, pairs = ret_down_sample_read_pairs(pairs, sample, down_samples[sample], outdir)
                
                proc_capture, proc = make_proc_capture(n)
                all_procs.append(proc)
                if write:
                    fh.write(cmd)
                    fh.write(proc_capture)
                
                n+=1
                lane_pairs[sample][lane] = pairs
    if write:
        wait_command = "wait %s\n" % " ".join(all_procs)
        fh.write(wait_command)
        
    return lane_pairs, scriptname
                
def ret_bwa_align_sort_cmd(pairs, sample, outdir, cpus, fasta_ref_path):
    cpus = str(cpus)
    r_pair = pairs[0].path
    f_pair = pairs[1].path
    lane_id = pairs[0].lane_id
    read_group = pairs[0].read_group

    output_bam = '%s/%s.bam' % (outdir, lane_id)
    bwa_command = "bwa mem -M -R \"%s\" -t %s %s %s %s" % (
        read_group, cpus, fasta_ref_path, r_pair, f_pair)
    sam_to_sorted_bam = ('samtools view -Sb -@ %s - | '
                         'samtools sort - -@ %s -o %s') % (
        cpus, cpus, output_bam)

    align_sort_cmd = "(%s | %s) &\n" % (bwa_command, sam_to_sorted_bam)

    return align_sort_cmd, output_bam


def ret_dedup_cmd(sample_input,  sample, outdir, picard_tools_path):
    dedup = '%s/MarkDuplicates.jar' % picard_tools_path
    dedup_basename = os.path.basename(sample_input).replace('.bam', '')
    dedup_metric = '%s/%s_metrics.txt' % (outdir, dedup_basename)
    dedup_output = '%s/%s.dedup.bam' % (outdir, dedup_basename)
    dedup_command = "java -Xmx4g -Djava.io.tmpdir=/scratch/dut -jar %s I=%s O=%s M=%s &\n" % (
        dedup, sample_input, dedup_output, dedup_metric)

    return dedup_command, dedup_output


def ret_merge_bams_cmd(sample_input,  sample, outdir, cpus, tempdir):
    cpus = str(cpus)
    bam_files = ' '.join(sample_input)
    temp_file = '%s/%s.temp' % (tempdir, sample)
    out_bam = '%s/%s.merged.bam' % (outdir, sample)
    merge_cmd = ('(samtools merge -@ %s - %s | samtools sort '
                 '- -m 10G -@ %s -T %s -o %s; samtools index %s) &\n') % (
                     cpus, bam_files, cpus, temp_file, out_bam, out_bam)

    return merge_cmd, out_bam


def ret_index_cmd(sample_input, sample, outdir):
    index_cmd = 'samtools index -b %s &\n' % sample_input

    return index_cmd, sample_input


def ret_realign_indels_cmd(sample_input, sample, outdir, cpus, fasta_ref_path, gatk_path):
    realign_basename = os.path.basename(sample_input).replace('.bam', '')
    realign_list = '%s/%s_target_intervals.list' % (outdir, realign_basename)
    out_bam = '%s/%s.realigned.bam' % (outdir, realign_basename)
    targets_cmd = '%s -T RealignerTargetCreator -R %s -I %s -o %s -nt %s' % (
        gatk_path, fasta_ref_path, sample_input, realign_list, cpus)

    raln_cmd = '%s -T IndelRealigner -R %s -I %s -targetIntervals %s -o %s' % (
        gatk_path, fasta_ref_path, sample_input, realign_list, out_bam)

    target_realign_cmd = '(%s; %s) &\n' % (targets_cmd, raln_cmd)

    return target_realign_cmd, out_bam


def ret_haplotype_caller_cmd(sample_input, sample, outdir, fasta_ref_path,
                             gatk_path, gatk_args, i, param_string=None):
    out_vars = '%s/%s.%s_raw_var.vcf' % (outdir, sample, str(i))
    if param_string:
        cmd = '%s -T HaplotypeCaller %s -R %s -I %s %s -o %s' % (
            gatk_path, param_string, fasta_ref_path,
            sample_input, gatk_args, out_vars)
    else:
        cmd = '%s -T HaplotypeCaller -R %s -I %s %s -o %s' % (
            gatk_path, fasta_ref_path, sample_input, gatk_args, out_vars)

    return cmd, out_vars


def ret_extract_snps_cmd(input_vars, sample, gatk_path, fasta_ref_path,
                         outdir, i):
    out_snps = '%s/%s.%s_raw_snps.vcf' % (outdir, sample, str(i))
    cmd = '%s -T SelectVariants -R %s -V %s -selectType SNP -o %s' % (
        gatk_path, fasta_ref_path, input_vars, out_snps)

    return cmd, out_snps


def ret_filter_snps_cmd(input_snps, sample, gatk_path, snp_filter,
                        snp_filter_name, fasta_ref_path, outdir, i):
    out_snps = '%s/%s.%s_filt_snps.vcf' % (outdir, sample, str(i))
    cmd = '%s -T VariantFiltration -R %s -V %s --filterExpression \"%s\" --filterName \"%s\" -o %s' % (
        gatk_path, fasta_ref_path, input_snps, snp_filter, snp_filter_name, out_snps)

    return cmd, out_snps


def ret_extract_indels_cmd(input_vars, sample, gatk_path, fasta_ref_path,
                           outdir, i):
    out_indels = '%s/%s.%s_raw_indels.vcf' % (outdir, sample, str(i))
    cmd = '%s -T SelectVariants -R %s -V %s -selectType INDEL -o %s' % (
        gatk_path, fasta_ref_path, input_vars, out_indels)

    return cmd, out_indels


def ret_filter_indels_cmd(input_indels, sample, gatk_path, indel_filter,
                          indel_filter_name, fasta_ref_path, outdir, i):
    out_indels = '%s/%s.%s_filt_indels.vcf' % (outdir, sample, str(i))
    cmd = '%s -T VariantFiltration -R %s -V %s --filterExpression \"%s\" --filterName \"%s\" -o %s' % (
        gatk_path, fasta_ref_path, input_indels, indel_filter,
        indel_filter_name, out_indels)

    return cmd, out_indels


def ret_combine_snps_cmd(all_snps, gatk_path, fasta_ref_path, outdir, i):
    out_snps = '%s/all_snps_%s.vcf' % (outdir, str(i))
    cmd = '%s -R %s -T CombineVariants --variant %s -o %s --excludeNonVariants --minimumN 1' % (
        gatk_path, fasta_ref_path, ' --variant '.join(all_snps), out_snps)

    return cmd, out_snps


def ret_combine_indels_cmd(all_indels, gatk_path, fasta_ref_path, outdir, i):
    out_indels = '%s/all_indels_%s.vcf' % (outdir, str(i))
    cmd = '%s -R %s -T CombineVariants --variant %s -o %s --excludeNonVariants --minimumN 1' % (
        gatk_path, fasta_ref_path, ' --variant '.join(all_indels), out_indels)

    return cmd, out_indels


def ret_vcf_tools_snp_cmd(input_snps, snp_filter_name):
    out_snps = input_snps.replace('.vcf', '')

    cmd = 'vcftools --vcf %s --remove-filtered LowQual --remove-filtered %s --recode --recode-INFO-all --out %s' % (
        input_snps, snp_filter_name, out_snps)

    out_snps = out_snps + '.recode.vcf'
    return cmd, out_snps


def ret_vcf_tools_indel_cmd(input_indels, indel_filter_name):
    out_indels = input_indels.replace('.vcf', '')
    cmd = 'vcftools --vcf %s --remove-filtered LowQual --remove-filtered %s --recode --recode-INFO-all --out %s' % (
        input_indels, indel_filter_name, out_indels)

    out_indels = out_indels + '.recode.vcf'
    return cmd, out_indels


def ret_first_recal_cmd(input_indels, input_snps, sample_input, sample, outdir,
                        gatk_path, fasta_ref_path, i):
    before_table = '%s/%s_%s.before_table' % (outdir, sample, str(i))
    cmd = '%s -T BaseRecalibrator -R %s -I %s -knownSites %s -knownSites %s -o %s' % (
        gatk_path, fasta_ref_path, sample_input, input_indels, input_snps,
        before_table)
    return cmd, before_table


def ret_print_reads_cmd(sample_input, recal_table, sample, outdir, gatk_path,
                        fasta_ref_path, i):
    out_bam = os.path.basename(sample_input).replace(
        '.bam', '.recal_%s.bam' % str(i))
    out_bam = '%s/%s' % (outdir, out_bam)
    cmd = '%s -T PrintReads -R %s -I %s -BQSR %s -o %s' % (
        gatk_path, fasta_ref_path, sample_input, recal_table, out_bam)

    return cmd, out_bam


def ret_second_recal_cmd(input_indels, input_snps, sample_input, before_table,
                         sample, outdir, gatk_path, fasta_ref_path, i):
    after_table = '%s/%s_%s.after_table' % (outdir, sample, str(i))
    cmd = '%s -T BaseRecalibrator -R %s -I %s -knownSites %s -knownSites %s -BQSR %s -o %s' % (
        gatk_path, fasta_ref_path, sample_input, input_snps, input_indels,
        before_table, after_table)
    return cmd, after_table


def ret_recal_plot_cmd(before_table, after_table, sample, outdir, gatk_path,
                       fasta_ref_path, i):
    recal_plot = '%s/%s_%s.plots' % (outdir, sample, str(i))
    cmd = '%s -T AnalyzeCovariates -R %s -before %s -after %s -plots %s' % (
        gatk_path, fasta_ref_path, before_table, after_table, recal_plot)

    return cmd


def ret_genotype_gvcf_command(sample_inputs, fasta_ref_path, gatk_path, cpus, outdir):
    all_inputs = '--variant %s' % ' --variant '.join(
        [item for sublist in sample_inputs.values() for item in sublist])
    output_gvcf = '%s/jg_%s.gvcf' % (outdir, '_'.join(sample_inputs.keys()))
    cmd = '%s -T GenotypeGVCFs -nt %s -ploidy 2 -R %s %s -allSites -o %s\n' % (
        gatk_path, cpus, fasta_ref_path, all_inputs, output_gvcf)

    return cmd, output_gvcf


def ret_hc_call_cmds(sample_input, sample, fasta_ref_path, gatk_path, outdir, i, param_string=None):
    gatk_args = " -stand_call_conf 30 -stand_emit_conf 30 -mbq 17 "
    snp_filter = 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0'
    indel_filter = 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0'
    indel_filter_name = 'default_indel_filter'
    snp_filter_name = 'default_snp_filter'
    haplotype_caller_cmd, out_vars = ret_haplotype_caller_cmd(
        sample_input=sample_input, sample=sample, outdir=outdir,
        fasta_ref_path=fasta_ref_path, gatk_path=gatk_path,
        gatk_args=gatk_args, i=i, param_string=param_string)

    extract_snps_cmd, raw_snps = ret_extract_snps_cmd(
        input_vars=out_vars, sample=sample, gatk_path=gatk_path,
        fasta_ref_path=fasta_ref_path, outdir=outdir, i=i)

    filter_snps_cmd, out_snps = ret_filter_snps_cmd(
        input_snps=raw_snps, sample=sample, gatk_path=gatk_path,
        snp_filter=snp_filter, snp_filter_name=snp_filter_name,
        fasta_ref_path=fasta_ref_path, outdir=outdir, i=i)

    extract_indels_cmd, raw_indels = ret_extract_indels_cmd(
        input_vars=out_vars, sample=sample, gatk_path=gatk_path,
        fasta_ref_path=fasta_ref_path, outdir=outdir, i=i)

    filter_indels_cmd, out_indels = ret_filter_indels_cmd(
        input_indels=raw_indels, sample=sample, gatk_path=gatk_path,
        indel_filter=indel_filter, indel_filter_name=indel_filter_name,
        fasta_ref_path=fasta_ref_path, outdir=outdir, i=i)

    hc_commands = '(%s) &\n' % '; '.join(
        [haplotype_caller_cmd, extract_snps_cmd, filter_snps_cmd,
         extract_indels_cmd, filter_indels_cmd])

    return hc_commands, out_snps, out_indels


def ret_recode_indels_and_snps_cmd(all_snps, all_indels, snp_filter_name,
                                   indel_filter_name, gatk_path,
                                   fasta_ref_path, outdir, i):

    combine_snps_cmd, out_snps = ret_combine_snps_cmd(
        all_snps=all_snps, gatk_path=gatk_path, fasta_ref_path=fasta_ref_path,
        outdir=outdir, i=i)

    combine_indels_cmd, out_indels = ret_combine_indels_cmd(
        all_indels=all_indels, gatk_path=gatk_path,
        fasta_ref_path=fasta_ref_path, outdir=outdir, i=i)

    vcf_tools_snp_cmd, recode_snps = ret_vcf_tools_snp_cmd(
        input_snps=out_snps, snp_filter_name=snp_filter_name)

    vcf_tools_indel_cmd, recode_indels = ret_vcf_tools_indel_cmd(
        input_indels=out_indels, indel_filter_name=indel_filter_name)

    combine_recode_cmd = '%s\n' % '\n'.join(
        [combine_snps_cmd, combine_indels_cmd, vcf_tools_snp_cmd,
         vcf_tools_indel_cmd])

    return combine_recode_cmd, recode_snps, recode_indels


def ret_recalibration_cmds(sample_input, sample, outdir, input_snps,
                           input_indels, gatk_path, fasta_ref_path, i):
    first_recal_cmd, before_table = ret_first_recal_cmd(
        input_indels=input_indels, input_snps=input_snps, sample_input=sample_input,
        sample=sample, outdir=outdir, gatk_path=gatk_path,
        fasta_ref_path=fasta_ref_path, i=i)

    print_reads_cmd, out_bam = ret_print_reads_cmd(
        sample_input=sample_input, recal_table=before_table, sample=sample,
        outdir=outdir, gatk_path=gatk_path, fasta_ref_path=fasta_ref_path,
        i=i)

    second_recal_cmd, after_table = ret_second_recal_cmd(
        input_indels=input_indels, input_snps=input_snps, sample_input=sample_input,
        before_table=before_table, sample=sample, outdir=outdir,
        gatk_path=gatk_path, fasta_ref_path=fasta_ref_path, i=i)

    recal_plot_cmd = ret_recal_plot_cmd(
        before_table=before_table, after_table=after_table, sample=sample,
        outdir=outdir, gatk_path=gatk_path, fasta_ref_path=fasta_ref_path, i=i)

    recal_cmd = '(%s) &\n' % '; '.join(
        [first_recal_cmd, print_reads_cmd, second_recal_cmd, recal_plot_cmd])

    return recal_cmd, out_bam


def final_variant_calling(input_dict, fasta_ref_path, gatk_path, cpus, write,
                          write_prefix, hcgvcfs):
    if write:
        fh, scriptname = make_script_fh(write_prefix)
        fh.write("#!/bin/bash\n")
    gatk_args = " -stand_call_conf 30 -stand_emit_conf 30 -mbq 17 "
    param_string = '--variant_index_type LINEAR --variant_index_parameter 128000 -ERC GVCF'
    hc_outdir = make_output_dir('final_variant_calling/hc_variant_calling')
    jg_outdir = make_output_dir('final_variant_calling/joint_genotypes')
    variant_dict = defaultdict(list)
    snp_filter = 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0'
    indel_filter = 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0'
    indel_filter_name = 'default_indel_filter'
    snp_filter_name = 'default_snp_filter'
    p = 0
    job_procs = []

    for sample in input_dict:
        cmd, out_var = ret_haplotype_caller_cmd(
            sample_input=input_dict[sample][0], sample=sample,
            outdir=hc_outdir, fasta_ref_path=fasta_ref_path,
            gatk_path=gatk_path, gatk_args=gatk_args,
            i='final', param_string=param_string)
        variant_dict[sample].append(out_var)
        proc_capture, proc = make_proc_capture(p)
        job_procs.append(proc)
        if write:
            fh.write(cmd + ' &\n')
            fh.write(proc_capture)
        p += 1

    wait_command = "wait %s\n" % " ".join(job_procs)
    job_procs = []
    
    if write:
        fh.write(wait_command)

    if hcgvcfs:
        for sample, gvcf in hcgvcfs:
            variant_dict[sample].append(gvcf)
        
    genotype_cmd, out_gvcf = ret_genotype_gvcf_command(
        sample_inputs=variant_dict, fasta_ref_path=fasta_ref_path,
        gatk_path=gatk_path, cpus=cpus, outdir=jg_outdir)

    if write:
        fh.write(genotype_cmd)
    sample = '_'.join(input_dict.keys())
    extract_snps_cmd, raw_snps = ret_extract_snps_cmd(
        input_vars=out_gvcf, sample=sample, gatk_path=gatk_path,
        fasta_ref_path=fasta_ref_path, outdir=jg_outdir, i='final')

    filter_snps_cmd, out_snps = ret_filter_snps_cmd(
        input_snps=raw_snps, sample=sample, gatk_path=gatk_path,
        snp_filter=snp_filter, snp_filter_name=snp_filter_name,
        fasta_ref_path=fasta_ref_path, outdir=jg_outdir, i='final')

    extract_indels_cmd, raw_indels = ret_extract_indels_cmd(
        input_vars=out_gvcf, sample=sample, gatk_path=gatk_path,
        fasta_ref_path=fasta_ref_path, outdir=jg_outdir, i='final')

    filter_indels_cmd, out_indels = ret_filter_indels_cmd(
        input_indels=raw_indels, sample=sample, gatk_path=gatk_path,
        indel_filter=indel_filter, indel_filter_name=indel_filter_name,
        fasta_ref_path=fasta_ref_path, outdir=jg_outdir, i='final')

    remove_nocall_cmd = 'cat %s | grep -v \'\./\.\'  > %s' % (
        out_gvcf, out_gvcf.replace('.gvcf', '.no_call_removed.gvcf'))
    badsnps_cmd = 'cat %s | grep \'%s\' | cut -f 1,2 > %s/bad_snps.tsv' % (
        out_snps, snp_filter_name, jg_outdir)
    badindels_cmd = 'cat %s   | grep \'%s\' | cut -f 1,2 > %s/bad_indels.tsv' % (
        out_indels, indel_filter_name, jg_outdir)
    bad_positions_cmd = 'cat %s/bad_indels.tsv %s/bad_snps.tsv > %s/bad_positions.tsv' % (
        jg_outdir, jg_outdir, jg_outdir)
    remove_bad_positions = 'vcftools --vcf %s --recode --out %s --exclude-positions %s/bad_positions.tsv' % (
        out_gvcf.replace('.gvcf', '.no_call_removed.gvcf'), out_gvcf.replace('.gvcf', '.no_call_removed.filt'), jg_outdir)

    if write:
        fh.write(extract_snps_cmd + '\n')
        fh.write(filter_snps_cmd + '\n')
        fh.write(extract_indels_cmd + '\n')
        fh.write(filter_indels_cmd + '\n')
        fh.write(remove_nocall_cmd + '\n')
        fh.write(badsnps_cmd + '\n')
        fh.write(badindels_cmd + '\n')
        fh.write(bad_positions_cmd + '\n')
        fh.write(remove_bad_positions + '\n')

        fh.close()

    return out_gvcf, scriptname


def training_commands(input_dict, rounds, fasta_ref_path, gatk_path, write,
                      write_prefix):
    variant_outdir = make_output_dir('training/variant_calling')
    recal_outdir = make_output_dir('training/recal_bams')
    if write:
        fh, scriptname = make_script_fh(write_prefix)
        fh.write('#!/bin/bash\n')
    job_procs = []
    p = 0
    print(rounds, range(rounds))
    for i in range(rounds):
        all_snps = []
        all_indels = []
        if i == 0:
            bam_dict = input_dict
        else:
            bam_dict = recal_bams
        for sample in bam_dict:
            hc_cmd, snps, indels = ret_hc_call_cmds(
                sample_input=bam_dict[sample][0], sample=sample,
                fasta_ref_path=fasta_ref_path, gatk_path=gatk_path,
                outdir=variant_outdir, i=i)

            all_snps.append(snps)
            all_indels.append(indels)
            proc_capture, proc = make_proc_capture(p)
            job_procs.append(proc)
            if write:
                fh.write(hc_cmd)
                fh.write(proc_capture)
            p += 1
        wait_command = "wait %s\n" % " ".join(job_procs)
        job_procs = []
        if write:
            fh.write(wait_command)

        recode_cmd, recode_snps, recode_indels = ret_recode_indels_and_snps_cmd(
            all_snps=all_snps, all_indels=all_indels,
            snp_filter_name='default_snp_filter',
            indel_filter_name='default_indel_filter', gatk_path=gatk_path,
            fasta_ref_path=fasta_ref_path, outdir=variant_outdir, i=i)

        if write:
            fh.write(recode_cmd)
        recal_bams = defaultdict(list)
        if i < rounds - 1:
            for sample in bam_dict:
                recal_cmd, recal_bam = ret_recalibration_cmds(
                    input_snps=recode_snps, input_indels=recode_indels,
                    sample_input=bam_dict[sample][0], sample=sample,
                    gatk_path=gatk_path, fasta_ref_path=fasta_ref_path,
                    outdir=recal_outdir, i=i)

                recal_bams[sample].append(recal_bam)
                proc_capture, proc = make_proc_capture(p)
                job_procs.append(proc)
                if write:
                    fh.write(recal_cmd)
                    fh.write(proc_capture)
                p += 1
            wait_command = "wait %s\n" % " ".join(job_procs)
            job_procs = []
            if write:
                fh.write(wait_command)
    fh.close()

    return recode_snps, recode_indels, scriptname


def parallel_cmd_loop(input_dict, nested, write, write_prefix,
                      cmd_function, outdir, **kwargs):
    out_dict = defaultdict(list)
    outdir = make_output_dir(outdir)
    if write:
        fh, scriptname = make_script_fh(write_prefix)
        fh.write("#!/bin/bash\n")
    n = 1
    all_procs = []
    job_procs = []
    out_dict = defaultdict(list)
    for sample in input_dict:
        if nested:
            try:
                nested_inputs = input_dict[sample].values()
            except AttributeError:
                nested_inputs = input_dict[sample]

            for nested_input in nested_inputs:
                cmd, output = cmd_function(
                    nested_input, sample, outdir, **kwargs)
                proc_capture, proc = make_proc_capture(n)
                if write:
                    fh.write(cmd)
                    fh.write(proc_capture)

                job_procs.append(proc)
                all_procs.append(proc)

                if n % 3 == 0:
                    wait_command = "wait %s\n" % " ".join(job_procs)
                    job_procs = []
                    if write:
                        fh.write(wait_command)

                out_dict[sample].append(output)
                # print(out_dict)
                n += 1
        else:
            # print(input_dict[sample])
            if type(input_dict[sample]) == list and len(input_dict[sample]) == 1:
                sample_input = input_dict[sample][0]
            else:
                sample_input = input_dict[sample]

            cmd, output = cmd_function(
                sample_input=sample_input, sample=sample, outdir=outdir, **kwargs)
            # print(output)
            proc_capture, proc = make_proc_capture(n)
            if write:
                fh.write(cmd)
                fh.write(proc_capture)

            job_procs.append(proc)
            all_procs.append(proc)

            if n % 3 == 0:
                wait_command = "wait %s\n" % " ".join(job_procs)
                job_procs = []
                if write:
                    fh.write(wait_command)

            out_dict[sample].append(output)

            n += 1

    if write:
        wait_command = "wait %s\n" % " ".join(all_procs)
        fh.write(wait_command)

    return out_dict, scriptname


def run(args):

    molng = lims.read_list_file(args.molngs)
    samples = lims.read_list_file(args.samples)
    flowcells = lims.read_list_file(args.flowcells)
    train = args.train
    
    fasta_path = args.fasta_path
    picard_tools_path = args.picard_path
    gatk = args.gatk_execute

    print(args.hcgvcfs)
    
    tigris_report = lims.Indexed_sample_reports(molng, samples, new_sample_report=args.new_sample_report)
    tigris_report.select_flowcells(flowcells)
    lane_pairs = tigris_report.ret_nested_pair_by_lane()
    scripts = []

    
    fasta_ref_path, s0 = prep_reference_commands(
        fasta_path, picard_tools_path, True)
    scripts.append(s0)

    if args.down_sample:
        print(args.down_sample)
        lane_pairs, dwn_s = down_sample_cmds(lane_pairs, args.down_sample, write=True,
                                            write_prefix='gatk_down_sample',outdir='preprocessing/down_sample')
        scripts.append(dwn_s)
        
    lane_bams, s1 = parallel_cmd_loop(
        input_dict=lane_pairs, nested=True, write=True,
        write_prefix='gatk_1_map_to_ref', cmd_function=ret_bwa_align_sort_cmd,
        outdir='preprocessing/map_to_ref_output', cpus=8,
        fasta_ref_path=fasta_ref_path)
    scripts.append(s1)
    
    dedup_lane_bams, s2 = parallel_cmd_loop(
        input_dict=lane_bams, nested=True, write=True,
        write_prefix='gatk_2_dedup_individual', cmd_function=ret_dedup_cmd,
        outdir='preprocessing/dedup_individuals',
        picard_tools_path=picard_tools_path)
    scripts.append(s2)
    
    sample_bams, s3 = parallel_cmd_loop(
        input_dict=dedup_lane_bams, nested=False, write=True,
        write_prefix='gatk_3_merge_bams', cmd_function=ret_merge_bams_cmd,
        outdir='preprocessing/merged_bams', cpus=8, tempdir='/scratch/dut')
    scripts.append(s3)
    
    deduped_sample_bams, s4 = parallel_cmd_loop(
        input_dict=sample_bams, nested=False, write=True,
        write_prefix='gatk_4_dedup_merged', cmd_function=ret_dedup_cmd,
        outdir='preprocessing/dedup_merged',
        picard_tools_path=picard_tools_path)
    scripts.append(s4)
    
    deduped_sample_bams, s5 = parallel_cmd_loop(
        input_dict=deduped_sample_bams, nested=False, write=True,
        write_prefix='gatk_5_index', cmd_function=ret_index_cmd,
        outdir='preprocessing/dedup_merged')
    scripts.append(s5)
    
    realigned_sample_bams, s6 = parallel_cmd_loop(
        input_dict=deduped_sample_bams, nested=False, write=True,
        write_prefix='gatk_6_realign_bams', cmd_function=ret_realign_indels_cmd,
        outdir='preprocessing/realigned_merged_bams', cpus=8,
        fasta_ref_path=fasta_ref_path, gatk_path=gatk)
    scripts.append(s6)
    
    if train:
        known_snps, known_indels, s7 = training_commands(
            input_dict=realigned_sample_bams, rounds=4, fasta_ref_path=fasta_ref_path,
            gatk_path=gatk, write=True, write_prefix='gatk_7_training')
        scripts.append(s7)
    elif args.known_snps and args.known_indels:
        known_snps = args.known_snps
        known_indels = args.known_indels
    else:
        print('train must be set to true or known snps and indels must be provieded')
        sys.exit(1)
        
    recalibrated_bams, s8 = parallel_cmd_loop(
        input_dict=realigned_sample_bams, nested=False, write=True,
        write_prefix='gatk_8_recalibrate_bams', outdir='final_variant_calling/recal_bams',
        cmd_function=ret_recalibration_cmds, input_snps=known_snps,
        input_indels=known_indels, gatk_path=gatk, fasta_ref_path=fasta_ref_path, i='final')
    scripts.append(s8)
    
    gvcf, s9 = final_variant_calling(input_dict=recalibrated_bams, fasta_ref_path=fasta_ref_path,
                                     gatk_path=gatk, cpus=24, write=True,
                                     write_prefix='gatk_9_final_variant_calling', hcgvcfs = args.hcgvcfs)
    scripts.append(s9)
    
    fh, _ = make_script_fh('master')

    fh.write('#!/bin/bash\n')
    
    commands = ['echo \"%s\"\ntime %s &> %s\n' % (
        x.replace('./', '').replace('.sh', ''),
        x,
        x.replace('./', '').replace('.sh', '.out'),
    ) for x in scripts]
    commands = ''.join(commands)
    fh.write(commands)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-molngs', action='store', dest='molngs',
                        help='path to file with all molngs paths seperated by newlines')

    parser.add_argument('-samples', action='store', dest='samples',
                        help='path to file with all sample names seperated by newlines')

    parser.add_argument('-flowcells', action='store', dest='flowcells',
                        help='path to file with all flowcell names seperated by newlines')

    parser.add_argument ('-down_sample', nargs=2, action='append',
                         help='pass sample name and fraction to down sample to, e.g -down_sample Atig_122 0.33')
    
    parser.add_argument('-cpus', action='store', dest='cpus', default='8',
                        help='number of cpus to use per sample')

    parser.add_argument('-fasta_path', action='store', dest='fasta_path',
                        default='/home/dut/projects/tigris/genome_annotation/fasta/tigris_scaffolds_filt_10000.fa',
                        help='path to genome fasta')
    
    parser.add_argument('-picard_path', action='store', dest='picard_path',
                        default='/home/dut/bin/picard-tools-1.119/',
                        help='path to picard tools')

    parser.add_argument('-gatk_execute', action='store', dest='gatk_execute',
                        default='java -Xmx3g -jar /home/dut/bin/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar',
                        help='gatk executable prefix')

    parser.add_argument('-train', action='store_true', dest='train',
                        default=False, help='train for known indels and snps')
    
    parser.add_argument('-known_snps', action='store', dest='known_snps',
                        help='path to known snps(required if training set to false)')

    parser.add_argument('-known_indels', action='store', dest='known_indels',
                        help='path to known indels(required if training set to false)')

    parser.add_argument('-add_hcgvcf',action='append',nargs=2, dest='hcgvcfs',
                        help='additional gvcf files produced by HaplotypeCaller to be used in joint genotypeing(-add_hcgvcf sample_name path_to_file')

    parser.add_argument('-new_sample_report', action='store_true', dest='new_sample_report',
                        default=False, help='train for known indels and snps')

    args = parser.parse_args()

    run(args)
