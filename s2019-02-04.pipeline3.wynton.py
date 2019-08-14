# Pipeline 3 imputes then does association testing for 1Mb chunks of simulated data

import numpy as np
import operator as op
import random
import math
import subprocess
import os
import pickle
import os.path
import gzip
import sys
from collections import defaultdict

# ===============================================
# INPUTS
# ===============================================
# vcf_filename_short = 'EUR-20000.5mb.vcf.gz'

# causal_bin_length = 10000
# num_causal_bins = 100

# region_start = 0
# region_end = 5000000

# ns = 20000
# nb = 1000000

# rho = 0.5
# tau = 0.5

# heritability = 0.4

# prevalence = 8 # top prevalence% of phenotypes are cases

# num_ref = 3000
# num_cc = 3000

# start_index = 0
# end_index = 1000000

# rvtest_set_length = 5000

# inputs from bash
region_start = int(sys.argv[1])
region_end = int(sys.argv[2])
ns = int(sys.argv[3])
nb = int(sys.argv[4])
causal_bin_length = int(sys.argv[5])
num_causal_bins = int(sys.argv[6])
rho = float(sys.argv[7])
tau = float(sys.argv[8])
heritability = float(sys.argv[9])
prevalence = float(sys.argv[10])
num_ref = int(sys.argv[11])
num_cc = int(sys.argv[12])
param_num = int(sys.argv[13])# this is necessary for base_name
temp_dir = str(sys.argv[14])
iteration_num = int(sys.argv[15])# is this really necessary? just use it in base_name, I guess?
p1_jobid = str(sys.argv[16])
population = str(sys.argv[17])
pop_size = str(sys.argv[18])
study_design = str(sys.argv[19])
start_index = int(sys.argv[20])
end_index = int(sys.argv[21])
rvtest_set_length = int(sys.argv[22])

vcf_filename_short = p1_jobid+'.'+population+'-'+pop_size+'.5mb.vcf.gz'

print('[INFO] Printing inputs to Python script read from bash')
print('region_start: '+str(region_start))
print('region_end :'+str(region_end))
print('ns: '+str(ns))
print('nb: '+str(nb))
print('causal_bin_length:'+str(causal_bin_length))
print('num_causal_bins: '+str(num_causal_bins))
print('rho: '+str(rho))
print('tau: '+str(tau))
print('heritability: '+str(heritability))
print('prevalence: '+str(prevalence))
print('num_ref: '+str(num_ref))
print('num_cc: '+str(num_cc))
print('param_num: '+str(param_num))
print('temp_dir: '+str(temp_dir))
print('iteration_num: '+str(iteration_num))
print('p1_jobid: '+str(p1_jobid))
print('population: '+str(population))
print('pop_size: '+str(pop_size))
print('study_design: '+str(study_design))
print('start_index: '+str(start_index))
print('end_index: '+str(end_index))
print('rvtest_set_length: '+str(rvtest_set_length))

data_path = temp_dir

print('[INFO] vcf_filename_short is: '+str(vcf_filename_short))

# ===============================================
# PATHS
# ===============================================
# mandrill_path = '/hernandez/mandrill/users/dmctong/imp_rvat/data'
mandrill_path = '/wynton/scratch/dmctong/imp_rvat/data'
support_path = mandrill_path+'/support'

# data_path = '/hernandez/mandrill/users/dmctong/imp_rvat/data/test' # temp_dir at scale
pheno_path = data_path+'/pheno'
rvtest_input_path = data_path+'/rvtest_input'
rvtest_output_path = data_path+'/rvtest_output'
result_path = data_path+'/result'

d_macs_file=data_path+'/'+vcf_filename_short[:-7]+'.causal_bin_length-'+str(causal_bin_length)+'.d_macs.txt' # moved from support_path 2018-09-17

subprocess.call('mkdir '+pheno_path, shell=True)
subprocess.call('mkdir '+rvtest_input_path, shell=True)
subprocess.call('mkdir '+rvtest_output_path, shell=True)
subprocess.call('mkdir '+result_path, shell=True)

# ===============================================
# FILE INPUTS + NAMING
# ===============================================
# base_name = 'test'
base_name = str(param_num)+'.'+str(iteration_num)
vcf_filename = mandrill_path+'/input-gzip/'+str(vcf_filename_short)
bgz_vcf_filename = mandrill_path+'/input-bgz/'+str(vcf_filename_short)

ad_filename = support_path+'/'+vcf_filename_short[:-7]+'.ad.gz' #additive model file name
vcf_mac_filename = support_path+'/'+vcf_filename_short[:-7]+'.mac.gz'
sfs_filename = support_path+'/'+vcf_filename_short[:-7]+'.d_sfs.gz'

array_sfs_file = support_path+'/InfiniumOmni2-5-8v1-3_A1_PopulationReport_MAFHistogram.txt'

mac_s_filename = support_path+'/'+str(population)+'.mac_s.all.sorted.gz'

# ===============================================
# IMPUTATION
# ===============================================

cc_ds_gen = rvtest_input_path+'/'+base_name+'.cc.ds.gen'
cc_ds_gen_sample = rvtest_input_path+'/'+base_name+'.cc.ds.gen.sample'

ref_haps = rvtest_input_path+'/'+base_name+'.ref.haps'
ref_legend = rvtest_input_path+'/'+base_name+'.ref.legend'
ref_sample = rvtest_input_path+'/'+base_name+'.ref.sample'

ref_haps_gz = ref_haps+'.gz'
ref_sample_gz = ref_sample+'.gz'
ref_legend_gz = ref_legend+'.gz'

# -----------------------------------------------------
# 1. Prephase using SHAPEIT
recomb_rate_file_name = '/netapp/home/dmctong/imp_rvat/scripts/genetic_map_chr22_combined_b37.20140701.txt'
phasing_output_file_name = rvtest_input_path+'/'+base_name+'.'+str(start_index)+'.cc.phasing.impute2'
Ne = 20000

subprocess.call('gzip '+cc_ds_gen,shell=True) # SHAPEIT works if you only zip the .gen file and not the .gen.sample file

gen_file_name = cc_ds_gen+'.gz'
gen_sample_file_name = cc_ds_gen_sample

# SHAPEIT
# ./shapeit -M $RECOMB_FILE -O $PHASING_OUT -T 4 --input-gen $GEN_FILE GEN_SAMPLE_FILE --effective-size 20000 --input-from $START --input-to $END
shapeit_call = '/netapp/home/dmctong/programs/shapeit/bin/shapeit '
dash_M = '-M '+recomb_rate_file_name+' '
dash_O = '-O '+phasing_output_file_name+' '
dash_T = '-T 4 '
dash_input = '--input-gen '+gen_file_name+' '+gen_sample_file_name+' '
dash_Ne = '--effective-size '+str(Ne)+' '
dash_start = '--input-from '+str(start_index)+' '
dash_end = '--input-to '+str(end_index)+' '
prephase_call = shapeit_call+dash_M+dash_O+dash_T+dash_input+dash_Ne+dash_start+dash_end

if os.path.isfile(phasing_output_file_name+'.haps'):
    print('[INFO] SHAPEIT prephase file exists. Skipping SHAPEIT prephase step.')
elif os.path.isfile(phasing_output_file_name+'.haps.gz'):
    print('[INFO] SHAPEIT prephase file (gzipped) exists. Skipping SHAPEIT prephase step.')
else:
    subprocess.call(prephase_call,shell=True)

# -----------------------------------------------------
# 2. Impute using IMPUTE4

imputed_gen_file_name = rvtest_input_path+'/'+base_name+'.'+str(start_index)+'.cc.imputed'

# IMPUTE4
# ./IMPUTE4 -h $REF_HAPS -l $REF_LEGEND -g $PHASING_OUT -m $RECOMB_FILE -int $START $END -no_maf_align -o $IMPUTED_OUT_GEN
# ./impute4 -h <file.hap.gz> -l <file.legend> -g <file> -m <file> -int <s> <e> -no_maf_align -o <file>
# note: need to check -g format (needs to be 1 row per SNP, first five columns ID1, ID2, position, allele1, allele2), two columns (haplotypes) per individual, 0,1 only

impute4_call = '/netapp/home/dmctong/programs/impute4/impute4.r265.2 '
dash_m = '-m '+recomb_rate_file_name+' '
dash_h = '-h '+ref_haps_gz+' '
dash_l = '-l '+ref_legend_gz+' '
dash_g = '-g '+phasing_output_file_name+'.haps ' #phased samples to be imputed
dash_int = '-int '+str(start_index)+' '+str(end_index)+' '
dash_Ne = '-Ne '+str(Ne)+' '
dash_o = '-o '+imputed_gen_file_name+' '
impute4_bash_command = impute4_call+dash_h+dash_l+dash_g+dash_m+dash_int+dash_o+'-no_maf_align'


if os.path.isfile(imputed_gen_file_name+'.gen'):
    print('[INFO] IMPUTE4 output gen file exists. Skipping IMPUTE4 step.')
elif os.path.isfile(imputed_gen_file_name+'.gen.gz'):
    print('[INFO] IMPUTE4 output gen file (gzipped) exists. Skipping IMPUTE4 step.')
else:
    subprocess.call(impute4_bash_command,shell=True)

# -----------------------------------------------------
# 3. Convert cc.imp.gen to cc.imp.vcf

def convert_gen_to_vcf(gen_file_name, gen_sample_file_name, output_vcf_file):
    bcftools_convert_call='/netapp/home/dmctong/programs/bcftools/bcftools convert '
    dash_gensample2vcf='--gensample2vcf '+gen_file_name+','+gen_sample_file_name+' '
    dash_o='-o '+output_vcf_file
    convert_gen_to_vcf_call=bcftools_convert_call+dash_gensample2vcf+dash_o
    subprocess.call(convert_gen_to_vcf_call,shell=True)

def bgzip_tabix(input_vcf):
    gzip_input_vcf=input_vcf+'.gz'
    bgzip_call='/netapp/home/dmctong/programs/htslib/bgzip -f '+input_vcf
    tabix_call='/netapp/home/dmctong/programs/htslib/tabix -fp vcf '+gzip_input_vcf
    subprocess.call(bgzip_call,shell=True)
    subprocess.call(tabix_call,shell=True)

# Convert gen to VCF and bgzip+tabix that VCF
imputed_vcf_file_name = imputed_gen_file_name+'.vcf'
if os.path.isfile(imputed_gen_file_name+'.vcf'):
    print('[INFO] Imputed VCF file already exists.')
else:
    print('[INFO] convert imputed.gen to .vcf')
    convert_gen_to_vcf(imputed_gen_file_name+'.gen', cc_ds_gen_sample, imputed_vcf_file_name) # bug with cc_ds_gen_sample?
    print('[INFO] bgzip_tabix imputed.vcf')
    bgzip_tabix(imputed_vcf_file_name)

# ===============================================
# RVAT
# ===============================================

# -----------------------------------------------------
# 1. Split genome into M bins. Will have 0 to M-1 bins (because Python)
print('[INFO] Split genome into M bins of length '+str(rvtest_set_length))

# list of start positions of each rvtest set bin
l_rvtest_set_start = range(region_start, region_end, rvtest_set_length)

# -----------------------------------------------------
# 2. Calculate start/end of each bin. Output this to rvtest_set_file

# format is # set1 \t 1:begin-end \n
rvtest_set_file = rvtest_input_path+'/'+base_name+'.rvtest_sets.txt'
outfile = open(rvtest_set_file,'w')

for start_pos in l_rvtest_set_start:
    set_num = start_pos/rvtest_set_length
    end_pos = start_pos + rvtest_set_length - 1
    out_str = 'set'+str(set_num)+'\t'+'1:'+str(start_pos)+'-'+str(end_pos)+'\n'
    outfile.write(out_str)

outfile.close()

# -----------------------------------------------------
# 3. Determine ve of each rvtest bin

l_rvtest_set_start = range(region_start, region_end, rvtest_set_length)

causal_sets_file = pheno_path+'/causal_sets_file.txt'
infile = open(causal_sets_file, 'r')

ve_causal_bin = 100./num_causal_bins

line = infile.readline().strip().split()
causal_start = int(line[1])
causal_end = int(line[2])

l_rvtest_ve=[0]*len(l_rvtest_set_start) # stores ve for each rvtest bin

for i, start_pos in enumerate(l_rvtest_set_start):
    if line == []:
        print('eof')
        break
    end_pos = start_pos + rvtest_set_length - 1
    print('start_pos: '+str(start_pos)+'\t'+'end_pos: '+str(end_pos)+'\t'+'causal_start: '+str(causal_start)+'\t'+'causal_end: '+str(causal_end))
    if end_pos <= causal_start:
        # next rvtest set
        l_rvtest_ve[i] = 0
        continue
    while start_pos >= causal_end:
        # next causal bin, and hold on this rvtest_set (so use a while loop)
        line = infile.readline().strip().split()
        if line == []:
            break
        causal_start = int(line[1])
        causal_end = int(line[2])
    if start_pos <= causal_start and end_pos > causal_start:
        # min(b-A, B-A)
        ve = min(end_pos - causal_start, causal_end - causal_start) * ve_causal_bin / (causal_end - causal_start)
        print('ve1: '+str(ve))
        l_rvtest_ve[i] = ve
    if start_pos > causal_start and start_pos < causal_end:
        # min(b-a, B-a)
        ve = min(end_pos - start_pos, causal_end - start_pos) * ve_causal_bin / (causal_end - causal_start)
        print('ve2: '+str(ve))
        l_rvtest_ve[i] = ve

# -----------------------------------------------------
# 4. Run rvtests (SKAT, SKAT-O, KBAC)
# 5. Compare rvtest results to p-value threshold, and output to standardized output file (same for loop)

def run_rvtests(input_tabix_vcf_file_name, rvtest_phenotype_file_name, rvtest_outfile_name, set_file_name, chromosome, start_index, end_index, test):
    rvtest_call='/netapp/home/dmctong/programs/rvtests/executable/rvtest '
    dash_inVcf='--inVcf '+input_tabix_vcf_file_name+' '
    dash_pheno='--pheno '+rvtest_phenotype_file_name+' '
    dash_out='--out '+rvtest_outfile_name+' '
    if test=='kbac':
        dash_kernel='--kernel kbac '
    elif test == 'skat':
        dash_kernel='--kernel skat '
    elif test == 'skato':
        dash_kernel='--kernel skato '
    elif test == 'firth':
        dash_kernel='--single firth '
    elif test == 'wald':
        dash_kernel='--single wald '
    dash_setFile='--setFile '+set_file_name+' '
    dash_rangeList='--rangeList '+str(chromosome)+':'+str(int(start_index))+'-'+str(int(end_index))
    kbac_bash_call=rvtest_call+dash_inVcf+dash_pheno+dash_out+dash_kernel+dash_setFile+dash_rangeList
    subprocess.call(kbac_bash_call,shell=True)

def output_standard_rvtest_result(infile_name, outfile_name, threshold, start_index, end_index, region_start, rvtest_set_size, l_r_ve, imp):
    correction=(start_index-region_start)/rvtest_set_size
    infile=open(infile_name,'r')
    outfile=open(outfile_name,'a')
    for i,line in enumerate(infile): # i counts from 0
        line=line.strip().split()
        if line[0]=='Range':
            continue
        set_num=line[0]
        set_info=line[1]
        set_start=line[1].split(':')[1].split('-')[0]
        set_end=line[1].split(':')[1].split('-')[1]
        num_var=line[2]
        num_poly_var=line[3]
        if test=='skat':
            if line[6]!='NA':
                p_value=float(line[6])
                if p_value<=threshold:
                    rvtest_result=1
                else:
                    rvtest_result=0
        elif test=='skato':
            if line[7]!='NA':
                p_value=float(line[7])
                if p_value<=threshold:
                    rvtest_result=1
                else:
                    rvtest_result=0
        elif test=='kbac':
            if line[5]!='NA':
                p_value=float(line[5])
                if p_value<=threshold:
                    rvtest_result=1
                else:
                    rvtest_result=0
        # print('i+correction: '+str(i+correction))
        # print('i')
        # print(i)
        # print('correction')
        # print(correction)
        # print('len(l_r_ve)')
        # print(len(l_r_ve))
        if imp == 1:
            ve = l_r_ve[i+correction-1] # -1 to compensate for header; correction to compensate for 1Mb-ness
        elif imp == 0: # if not imputed, don't need correction since rvtests analyzes whole 5Mb region anyway (why???)
            ve = l_r_ve[i-1]
        print([test,imp,set_start,set_end,num_var,num_poly_var,p_value,ve,rvtest_result])
        if int(set_start) < int(end_index):
            outfile.write(str(test)+'\t'+str(imp)+'\t'+str(set_start)+'\t'+str(set_end)+'\t'+str(num_var)+'\t'+str(num_poly_var)+'\t'+str(p_value)+'\t'+str(ve)+'\t'+str(rvtest_result)+'\n')
        else:
            break
    infile.close()
    outfile.close()
    return 0

threshold = 0.05 / len(l_rvtest_set_start)

l_tests=['skat','kbac','skato'] # list of tests to run, don't forget SKAT-O later (for rvtests: skato)

# existing previous files
rvtest_pheno_filename = pheno_path+'/'+str(base_name)+'.rvtests.pheno.gz'

# imputed files
rvtest_imputed_out_file = rvtest_output_path+'/'+base_name+'.'+str(int(start_index))+'.cc.imp'
imputed_vcf_file_name_gz = rvtest_input_path+'/'+base_name+'.'+str(int(start_index))+'.cc.imputed.vcf.gz'

# print('imputed_vcf_file_name_gz: '+str(imputed_vcf_file_name_gz))

# full files
rvtest_full_out_file = rvtest_output_path+'/'+base_name+'.'+str(int(start_index))+'.cc.full'
cc_vcf_gz = rvtest_input_path+'/'+base_name+'.cc.vcf.gz'

rvtest_standard_out = rvtest_output_path+'/'+base_name+'.'+str(int(start_index))+'.rvtest.data.all'
subprocess.call('touch '+rvtest_standard_out,shell=True)

for test in l_tests:
    print(test)
    if test=='skat': # map test to the output file name (super annoying!)
        cap_test='Skat'
    elif test=='skato':
        cap_test='SkatO'
    elif test == 'kbac':
        cap_test='Kbac'
    # imputed data
    rvtest_imputed_out_test = rvtest_imputed_out_file+'.'+str(cap_test)+'.assoc'
    if not os.path.isfile(rvtest_imputed_out_test):
        run_rvtests(imputed_vcf_file_name_gz, rvtest_pheno_filename, rvtest_imputed_out_file, rvtest_set_file, 1, start_index, end_index, test)
    else:
        print('rvtest_imputed_out_test already exists!')
    rvtest_imputed_out_test = rvtest_imputed_out_file+'.'+str(cap_test)+'.assoc'
    output_standard_rvtest_result(rvtest_imputed_out_test, rvtest_standard_out, threshold, start_index, end_index, region_start, rvtest_set_length, l_rvtest_ve, imp=1)
    # full data
    rvtest_full_out_test = rvtest_full_out_file+'.'+str(cap_test)+'.assoc'
    if not os.path.isfile(rvtest_full_out_test):
        run_rvtests(cc_vcf_gz, rvtest_pheno_filename, rvtest_full_out_file, rvtest_set_file, 1, start_index, end_index, test)
    else:
        print('rvtest_full_out_test already exists!')
    rvtest_full_out_test = rvtest_full_out_file+'.'+str(cap_test)+'.assoc'
    output_standard_rvtest_result(rvtest_full_out_test, rvtest_standard_out, threshold, start_index, end_index, region_start, rvtest_set_length, l_rvtest_ve, imp=0)

# -----------------------------------------------------
# 6. Run Firth bias-corrected single variant association test

test = 'firth'

rvtest_imputed_out_file = rvtest_output_path+'/'+base_name+'.'+str(int(start_index))+'.cc.imp'
rvtest_full_out_file = rvtest_output_path+'/'+base_name+'.'+str(int(start_index))+'.cc.full'

# imputed
run_rvtests(imputed_vcf_file_name_gz, rvtest_pheno_filename, rvtest_imputed_out_file, rvtest_set_file, 1, start_index, end_index, test)
# full
run_rvtests(cc_vcf_gz, rvtest_pheno_filename, rvtest_full_out_file, rvtest_set_file, 1, region_start, region_end, test)
# this just seems to test them all, so screw the start_index part of the file name

# -----------------------------------------------------
# 7. Calculate set-based Firth output

# INPUT
# firth_(full or imp)_out       Range/Chrom/Pos/Ref/Alt/N_informative/Test/Beta/SE/P-value

# OUTPUT
# firth_set_out                 test/imp/set_start/set_end/num_vars/lowest p_value/ve/Firth_result

# input files
firth_full_out = rvtest_output_path+'/'+base_name+'.'+str(int(start_index))+'.cc.full.SingleFirth.assoc'
firth_imp_out = rvtest_output_path+'/'+base_name+'.'+str(int(start_index))+'.cc.imp.SingleFirth.assoc'

# output files
firth_set_out = rvtest_output_path+'/'+base_name+'.'+str(int(start_index))+'.cc.firth.data.all'

outfile = open(firth_set_out, 'w')

# full
test = 'firth'
imp = 0

infile = open(firth_full_out, 'r')
header = infile.readline()

v_num_vars = np.zeros(len(l_rvtest_set_start))
v_lowest_p = np.ones(len(l_rvtest_set_start))
v_firth_result = np.zeros(len(l_rvtest_set_start))

flag = True
for line in infile:
    line = line.strip().split()
    if line[9]=='NA':
        continue
    else:
        p_value = float(line[9])
    pos = int(line[2])
    if pos < start_index or pos > end_index: # we only want the 1Mb chunk right now
        break
    # rvtest_set = pos / rvtest_set_length # should be rvtest_set = (pos - region_start)/rvtest_set_length?
    rvtest_set = (pos - region_start)/rvtest_set_length
    v_num_vars[rvtest_set] = v_num_vars[rvtest_set] + 1
    if p_value <= v_lowest_p[rvtest_set]:
        v_lowest_p[rvtest_set] = p_value
    if p_value <= 5e-8:
        v_firth_result[rvtest_set] = 1

for i in range((start_index - region_start) / rvtest_set_length, (end_index - region_start) / rvtest_set_length):
    print('firth', str(imp), l_rvtest_set_start[i], l_rvtest_set_start[i]+rvtest_set_length-1, v_num_vars[i], v_lowest_p[i], l_rvtest_ve[i], v_firth_result[i])
    outfile.write(str(test)+'\t'+str(imp)+'\t'+str(l_rvtest_set_start[i])+'\t'+str(l_rvtest_set_start[i]+rvtest_set_length-1)+'\t'+str(v_num_vars[i])+'\t'+str(v_num_vars[i])+'\t'+str(v_lowest_p[i])+'\t'+str(l_rvtest_ve[i])+'\t'+str(v_firth_result[i])+'\n')

# imputed
test = 'firth'
imp = 1

infile = open(firth_imp_out, 'r')
header = infile.readline()

v_num_vars = np.zeros(len(l_rvtest_set_start))
v_lowest_p = np.ones(len(l_rvtest_set_start))
v_firth_result = np.zeros(len(l_rvtest_set_start))

flag = True
for line in infile:
    line = line.strip().split()
    if line[9]=='NA':
        continue
    else:
        p_value = float(line[9])
    pos = int(line[2])
    # rvtest_set = pos / rvtest_set_length
    rvtest_set = (pos - region_start)/rvtest_set_length
    v_num_vars[rvtest_set] = v_num_vars[rvtest_set] + 1
    if p_value <= v_lowest_p[rvtest_set]:
        v_lowest_p[rvtest_set] = p_value
    if p_value <= 5e-8:
        v_firth_result[rvtest_set] = 1

for i in range((start_index - region_start) / rvtest_set_length, (end_index - region_start) / rvtest_set_length):
    outfile.write(str(test)+'\t'+str(imp)+'\t'+str(l_rvtest_set_start[i])+'\t'+str(l_rvtest_set_start[i]+rvtest_set_length-1)+'\t'+str(v_num_vars[i])+'\t'+str(v_num_vars[i])+'\t'+str(v_lowest_p[i])+'\t'+str(l_rvtest_ve[i])+'\t'+str(v_firth_result[i])+'\n')

outfile.close()

# Delete ${PARAM_NUM}.${ITERATION_NUM}.${START_INDEX}.cc.imputed.gen
rm_call = 'rm '+imputed_gen_file_name+'.gen'
subprocess.call(rm_call, shell=True)


# ===============================================
# GWAS
# ===============================================

# Unzip rvtest.pheno file (it isn't too big)
subprocess.call('gunzip '+str(rvtest_pheno_filename),shell=True)

# Recode the rvtests.pheno file to use 0,1 instead of 1,2
rvtest_pheno_filename_unzip = pheno_path+'/'+str(base_name)+'.rvtests.pheno'
rvtest_pheno_recode_filename = rvtest_pheno_filename+'.recode'

infile = open(rvtest_pheno_filename_unzip,'r')
outfile = open(rvtest_pheno_recode_filename,'w')
infile.readline()
outfile.write('FID\tIID\tfatid\tmatid\tsex\ty1\ty2')
for line in infile:
    line = line.strip().split()
    if line[5]=='1':
        new_line=[line[0],line[1],'1']
    if line[5]=='2':
        new_line=[line[0],line[1],'2']
    outfile.write('\t'.join(new_line)+'\n')

outfile.close()
infile.close()

# Rezip rvtest.pheno file to save space
subprocess.call('gzip '+str(rvtest_pheno_filename_unzip),shell=True)

# Run PLINK GWAS on imputed files
# what is vcf? imputed vcf file (remember P3 only runs on 1Mb of 5Mb chunk)
# imputed_vcf_file_name_gz already exists from before

# what is out?
# let's make a GWAS output directory?
# so like gwas_output_directory+'/'+str(base_name)+'.'+str(start_index)

gwas_path = data_path+'/gwas'
subprocess.call('mkdir '+gwas_path, shell=True)

imputed_gwas_output_prefix = gwas_path+'/'+str(base_name)+'.'+str(start_index)+'.cc.imputed'

plink_call = '/netapp/home/dmctong/programs/plink/plink --vcf '+str(imputed_vcf_file_name_gz)+' --pheno '+str(rvtest_pheno_recode_filename)+' --assoc --allow-no-sex --double-id --maf 0.01 --out '+str(imputed_gwas_output_prefix)

subprocess.call(plink_call, shell=True)

# PLINK call looks like
# ./plink --vcf /hernandez/mandrill/users/dmctong/test/1.1.16000000.cc.imputed.vcf.gz --pheno /hernandez/mandrill/users/dmctong/test/1.1.rvtests.pheno --assoc --allow-no-sex --double-id --out /hernandez/mandrill/users/dmctong/test/stuff --maf 0.01

# Run PLINK GWAS on full files
# so full VCF file and change the output prefix

cc_vcf_gz = rvtest_input_path+'/'+base_name+'.cc.vcf.gz' # is full VCF file
full_gwas_output_prefix = gwas_path+'/'+str(base_name)+'.cc.full'

plink_call = '/netapp/home/dmctong/programs/plink/plink --vcf '+str(cc_vcf_gz)+' --pheno '+str(rvtest_pheno_recode_filename)+' --assoc --allow-no-sex --double-id --maf 0.01 --out '+str(full_gwas_output_prefix)
subprocess.call(plink_call, shell=True)

# Zip rvtest.pheno.recode file to save space
subprocess.call('gzip '+str(rvtest_pheno_recode_filename),shell=True)

# Process GWAS output files (.assoc) into existing standard format and store
# NB we will concatenate all standard formats into one file for plotting later
# output file format:
# TEST IMP SET_BEGIN SET_END NUM_VARS OR LOWEST_P VE RESULT

# imputed
test = 'gwas'
imp = 1

gwas_imp_out = imputed_gwas_output_prefix+'.assoc'

infile = open(gwas_imp_out, 'r')
header = infile.readline()

# output files
gwas_set_out = gwas_path+'/'+base_name+'.'+str(int(start_index))+'.cc.gwas.data.all'
outfile = open(gwas_set_out, 'w')

v_num_vars = np.zeros(len(l_rvtest_set_start))
v_lowest_p = np.ones(len(l_rvtest_set_start))
v_or = np.ones(len(l_rvtest_set_start))
v_firth_result = np.zeros(len(l_rvtest_set_start))

flag = True
for line in infile:
    line = line.strip().split()
    if line[0]=='CHR':
        continue
    if line[8]=='NA':
        continue
    else:
        p_value = float(line[8])
        if line[9] != 'NA':
            odds_ratio = float(line[9])
        else:
            odds_ratio = 999. #we don't really care about OR, so just put something in there
    pos = int(line[2])
    # rvtest_set = pos / rvtest_set_length
    rvtest_set = (pos - region_start)/rvtest_set_length
    v_num_vars[rvtest_set] = v_num_vars[rvtest_set] + 1
    if p_value <= v_lowest_p[rvtest_set]:
        v_lowest_p[rvtest_set] = p_value
        v_or[rvtest_set] = odds_ratio
    if p_value <= 5e-8:
        v_firth_result[rvtest_set] = 1

for i in range((start_index - region_start) / rvtest_set_length, (end_index - region_start) / rvtest_set_length):
    outfile.write(str(test)+'\t'+str(imp)+'\t'+str(l_rvtest_set_start[i])+'\t'+str(l_rvtest_set_start[i]+rvtest_set_length-1)+'\t'+str(v_num_vars[i])+'\t'+str(v_or[i])+'\t'+str(v_lowest_p[i])+'\t'+str(l_rvtest_ve[i])+'\t'+str(v_firth_result[i])+'\n')

# outfile.close()

# full
test = 'gwas'
imp = 0

gwas_full_out = full_gwas_output_prefix+'.assoc'

infile = open(gwas_full_out, 'r')
header = infile.readline()

v_num_vars = np.zeros(len(l_rvtest_set_start))
v_lowest_p = np.ones(len(l_rvtest_set_start))
v_or = np.ones(len(l_rvtest_set_start))
v_firth_result = np.zeros(len(l_rvtest_set_start))

flag = True
for line in infile:
    line = line.strip().split()
    if line[0]=='CHR':
        continue
    if line[8]=='NA': # p_value
        continue
    if line[9] == 'NA': # odds_ratio - why is this NA? not sure
        continue
    else:
        p_value = float(line[8])
        if line[9] != 'NA':
            odds_ratio = float(line[9])
        else:
            odds_ratio = 999. #we don't really care about OR, so just put something in there
    pos = int(line[2])
    if pos < start_index or pos > end_index: # we only want the 1Mb chunk right now
        continue
    # rvtest_set = pos / rvtest_set_length # should be rvtest_set = (pos - region_start)/rvtest_set_length?
    rvtest_set = (pos - region_start)/rvtest_set_length
    v_num_vars[rvtest_set] = v_num_vars[rvtest_set] + 1
    if p_value <= v_lowest_p[rvtest_set]:
        v_lowest_p[rvtest_set] = p_value
        v_or[rvtest_set] = odds_ratio
    if p_value <= 5e-8:
        v_firth_result[rvtest_set] = 1

for i in range((start_index - region_start) / rvtest_set_length, (end_index - region_start) / rvtest_set_length):
    print('gwas', str(imp), l_rvtest_set_start[i], l_rvtest_set_start[i]+rvtest_set_length-1, v_num_vars[i], v_lowest_p[i], l_rvtest_ve[i], v_firth_result[i])
    outfile.write(str(test)+'\t'+str(imp)+'\t'+str(l_rvtest_set_start[i])+'\t'+str(l_rvtest_set_start[i]+rvtest_set_length-1)+'\t'+str(v_num_vars[i])+'\t'+str(v_or[i])+'\t'+str(v_lowest_p[i])+'\t'+str(l_rvtest_ve[i])+'\t'+str(v_firth_result[i])+'\n')

outfile.close()


# DELETE unnecessary files
# ref.haps.gz
subprocess.call('rm '+ref_haps_gz,shell=True)
# ref.legend.gz
subprocess.call('rm '+ref_legend_gz,shell=True)
# ref.sample.gz
subprocess.call('rm '+ref_sample_gz,shell=True)
# cc.imputed.vcf.gz
subprocess.call('rm '+imputed_vcf_file_name_gz,shell=True)
# cc.phasing.impute2.haps
subprocess.call('rm '+phasing_output_file_name+'.haps',shell=True)
# cc.phasing.impute2.sample
subprocess.call('rm '+phasing_output_file_name+'.sample',shell=True)
