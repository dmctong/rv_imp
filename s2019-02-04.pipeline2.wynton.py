# Pipeline 2 simulates phenotypes for all N individuals

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

# Testing with iqint parameters

# vcf_filename_short = 'EUR-20000.5mb.vcf.gz'
# region_start = 16000000
# region_end = 21000000
# ns = 20000
# nb = 1000000
# causal_bin_length = 10000
# num_causal_bins = 10
# rho = 0.5
# tau = 0.5
# heritability = 0.01
# prevalence = 8 # top prevalence% of phenotypes are cases
# num_ref = 3000
# num_cc = 3000
# param_num = 1
# temp_dir = '/hernandez/mandrill/users/dmctong/imp_rvat/data/test' # temp_dir at scale
# param_num = 1
# iteration_num = 1

# inputs from bash
# vcf_filename_short = 'EUR-20000.5mb.vcf.gz'
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

vcf_filename_short = p1_jobid+'.'+population+'-'+pop_size+'.5mb.vcf.gz'

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
print('iteration_num: '+str(iteration_num))

data_path = temp_dir

# ===============================================
# PATHS
# ===============================================
# mandrill_path = '/hernandez/mandrill/users/dmctong/imp_rvat/data'
wynton_path = '/wynton/scratch/dmctong/imp_rvat/data'
support_path = wynton_path+'/support'

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
vcf_filename = wynton_path+'/input-gzip/'+str(vcf_filename_short) #TODO rsync from here to temp_dir?
bgz_vcf_filename = wynton_path+'/input-bgz/'+str(vcf_filename_short) #TODO rsync from here to temp_dir?

if iteration_num%2==0: # distribute across 2 files to reduce read/write errors
    ad_filename = support_path+'/'+vcf_filename_short[:-7]+'.ad.even.gz' #additive model file name
else:
    ad_filename = support_path+'/'+vcf_filename_short[:-7]+'.ad.odd.gz' #additive model file name
vcf_mac_filename = support_path+'/'+vcf_filename_short[:-7]+'.mac.gz'
sfs_filename = support_path+'/'+vcf_filename_short[:-7]+'.d_sfs.gz'

array_sfs_file = support_path+'/InfiniumOmni2-5-8v1-3_A1_PopulationReport_MAFHistogram.txt'

mac_s_filename = support_path+'/'+str(population)+'.mac_s.all.sorted.gz'

# ===============================================
# MAIN
# ===============================================

# -----------------------------------------------------
# 1. Split genome into M bins. Will have 0 to M-1 bins (because Python)
print('[INFO] Split genome into M bins of length '+str(causal_bin_length))

region_length = region_end - region_start

# list of start positions of each potential causal bin
l_pot_causal_start = range(region_start, region_end, causal_bin_length)

# -----------------------------------------------------
# 2. Pick K causal bins from M bins.
print('[INFO] Picking '+str(num_causal_bins)+'from M bins of length '+str(causal_bin_length))
l_causal_start = random.sample(l_pot_causal_start, num_causal_bins) # list of causal bin starting positions

# -----------------------------------------------------
# 3. Calculate start/end of each bin. Output this to causal_sets_file

causal_sets_file = pheno_path+'/causal_sets_file.txt'
outfile = open(causal_sets_file, 'w')

l_causal_start.sort()

for start_pos in l_causal_start:
    causal_bin_number=start_pos/causal_bin_length
    out_str = '1'+'\t'+str(start_pos)+'\t'+str(start_pos+causal_bin_length-1)+'\t'+str(causal_bin_number)+'\n'
    outfile.write(out_str)

outfile.close()

# -----------------------------------------------------
# 4. Create additive model file and (position,MAC) file (do once per VCF file)

# NB filenames are defined at the top of this file

if not os.path.isfile(ad_filename):
    if not os.path.isfile(vcf_mac_filename):
        print('[INFO] ad_filename does not exist. vcf_mac_filename does not exist. Generating ad_filename and vcf_mac_filename now.')
        outfile_ad = gzip.open(ad_filename, 'wb')
        outfile_mac = gzip.open(vcf_mac_filename, 'wb')
        # open the VCF file
        if vcf_filename.split('/')[-1].split('.')[-1]=='gz':
            vcf_file = gzip.open(vcf_filename, 'rb')
        # go through VCF file and find additive model
        for j,line in enumerate(vcf_file):
            if j % 10000 == 0: #207451 sites
                print('[INFO] additive model vcf file line: '+str(j))
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                line = line.strip().split()
                l_ids = line[9:len(line)]
            else:
                line = line.strip().split()
                position = line[1]
                # for each line, calculate the additive model and write to file
                l_out = []
                mac = 0
                for i, indiv in enumerate(line):
                    if i <= 8: # skip the id, ref, alt type of things
                        continue
                    else:
                        ad = [int(a) for a in indiv.split('|')]
                        mac = mac + sum(ad)
                        l_out.append(str(sum(ad)))
                outfile_ad.write(str(position)+'\t'+str(mac)+'\t'+'\t'.join(l_out)+'\n')
                outfile_mac.write(str(position)+'\t'+str(mac)+'\n')
        # close outfiles
        outfile_ad.close()
        outfile_mac.close()

if os.path.isfile(ad_filename): # if ad_filename already exists, we need the list of ID's (l_ids)
    if os.path.isfile(vcf_mac_filename):
        print('[INFO] ad_filename and vcf_mac_filename already exist')
        print('[INFO] finding l_ids in vcf_file now')
        # open the VCF file
        if vcf_filename.split('/')[-1].split('.')[-1]=='gz':
            vcf_file = gzip.open(vcf_filename, 'rb')
        # find the l_ids in the VCF file then break
        for j,line in enumerate(vcf_file):
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                line = line.strip().split()
                l_ids = line[9:len(line)]
            else:
                break

d_sfs = defaultdict(list) #d_sfs[MAC]=[list of positions]

if not os.path.isfile(sfs_filename):
    print('[INFO] sfs_filename does not yet exist. Generating d_sfs and outputting to sfs_filename now.')
    infile_mac = gzip.open(vcf_mac_filename, 'rb')
    for line in infile_mac:
        line = line.strip().split()
        mac = int(line[1])
        position = int(line[0])
        d_sfs[mac].append(position)
    outfile_sfs = gzip.open(sfs_filename, 'wb')
    # output sfs
    for mac in d_sfs:
        outfile_sfs.write(str(mac)+'\t'+'\t'.join([str(a) for a in d_sfs[mac]])+'\n')
    outfile_sfs.close()
else: # read in d_sfs if it already exists
    print('[INFO] sfs_filename already exists. Reading in d_sfs now.')
    infile_sfs = gzip.open(sfs_filename, 'rb')
    for line in infile_sfs:
        line = line.strip().split()
        mac = int(line[0])
        d_sfs[mac]=[int(a) for a in line[1:]]

# -----------------------------------------------------
# 5. Get (MAC,s) for every SNP in genome using SFS_CODE outfile. 
#     DATA: d_macs[MAC]=[list of s]
#     *not every MAC will be filled

# 5a.   If pickled d_macs does not exist, then get d_macs from SFS_CODE outfile.
#       If pickled d_macs does exist, load into d_macs

def read_acs(sfs_coder_file, ns, nb):
    l_sfs=[0]*(2*ns)
    d_macs=defaultdict(list)
    # if sfs_coder_file.split('.')[-1]=='gz':
    #     mac_s_file=gzip.open(sfs_coder_file,'rb')
    # if sfs_coder_file.split('.')[-1]=='vcf':
    mac_s_file=gzip.open(sfs_coder_file,'rb')
    for line in mac_s_file:
        line=line.strip().split()
        mac=int(line[0])
        sel_coeff=float(line[1])
        if sel_coeff==0 or mac==0:
            continue
        new_mac=np.random.binomial(2*ns,float(mac)/(2*nb),size=1)[0]
        if new_mac<2*ns and new_mac>0:
            d_macs[new_mac].append(sel_coeff)
    return d_macs

d_macs = defaultdict(list) # d_macs[mac] = [list of s]
if os.path.isfile(d_macs_file):
    print('[INFO] loading d_macs')
    d_macs = pickle.load( open(d_macs_file, "rb" ) ) # read from file
else:
    print('[INFO] reading d_macs from mac_s_filename')
    d_macs = read_acs(mac_s_filename, ns, nb)
    pickle.dump(d_macs, open(d_macs_file,'wb'))

# 5b.   Get (MAC, s) pair for each position in genome and output to d_macs_actual[mac] = [list of s]
from bisect import bisect_left
def takeClosest(myList, myNumber):
    """
    Assumes myList is sorted. Returns closest value to myNumber.

    If two numbers are equally close, return the smallest number.
    """
    pos = bisect_left(myList, myNumber)
    if pos == 0:
        return myList[0]
    if pos == len(myList):
        return myList[-1]
    before = myList[pos - 1]
    after = myList[pos]
    if after - myNumber < myNumber - before:
       return after
    else:
       return before

d_macs_actual = defaultdict(list) # d_macs_actual[mac] = [list of s] in simulated population
l_macs = d_macs.keys()
l_macs.sort()

print('[INFO] getting (MAC,s) pairs for all SNPs in input population')
infile = gzip.open(vcf_mac_filename, 'rb')
for i, line in enumerate(infile):
    if i % 1000 == 0:
        print('line: '+str(i))
    line = line.strip().split()
    position = int(line[0])
    mac = int(line[1])
    if mac in l_macs:
        s = random.sample(d_macs[mac], 1)[0]
    else:
        s = random.sample(d_macs[takeClosest(l_macs, mac)], 1)[0]
    d_macs_actual[mac].append(s)

with open(data_path+'/d_macs_actual.txt', 'wb') as handle: # what is this? where should it actually go?
    pickle.dump(d_macs_actual, handle, protocol=pickle.HIGHEST_PROTOCOL)

# -----------------------------------------------------
# 6. Determine effect sizes for each SNP within causal bins
print('[INFO] Getting effect sizes for each causal SNP')

N = len(l_ids)
infile = gzip.open(ad_filename, 'rb')
l_macs_actual = d_macs_actual.keys()
l_macs_actual.sort()
v_genetic_effect = np.zeros(N)

d_causals = defaultdict(float)

causal_start_counter = 0

l_flat_s = [y for x in d_macs_actual.values() for y in x] # d_macs_actual[mac] = [list of s] in simulated population; l_flat_s is just a list of all selection coefficients

def get_effect_size(position, mac, tau, rho, d_macs_actual, l_macs_actual, v_genetic_effect, v_snp, d_causals, l_flat_s):
    # pick delta
    delta = np.random.choice([-1,1],replace=False)
    # pick s
    if np.random.uniform()<=rho:
        if mac in l_macs_actual:
            s = np.random.choice(d_macs_actual[mac])
        else:
            s = np.random.choice(d_macs_actual[takeClosest(l_macs_actual, mac)])
        z = delta  * (abs(s) ** tau)
    else:
        s_r = random.choice(l_flat_s)
        z = delta * (math.fabs(s_r) ** tau)
    v_genetic_effect = v_genetic_effect + np.multiply(z, v_snp)
    d_causals[position] = z
    return [v_genetic_effect, d_causals]

for i, line in enumerate(infile): # infile is ad_filename
    if i % 1000 == 0:
        print(i)
    line = line.strip().split()
    position = int(line[0])
    mac = int(line[1])
    v_snp = np.asarray([int(a) for a in line[2:]])
    if (position >= l_causal_start[causal_start_counter] and position < l_causal_start[causal_start_counter]+causal_bin_length):
        [v_genetic_effect, d_causals] = get_effect_size(position, mac, tau, rho, d_macs_actual, l_macs_actual, v_genetic_effect, v_snp, d_causals, l_flat_s)
    elif position > l_causal_start[causal_start_counter] + causal_bin_length:
        causal_start_counter = causal_start_counter + 1
        print('causal_start_counter: '+str(causal_start_counter))
        if causal_start_counter > len(l_causal_start)-1:
            print('[WARN] Last causal bin reached. Breaking in effect size loop')
            break
        print('l_causal_start position: '+str(l_causal_start[causal_start_counter]))
        print('position: '+str(position))
        if (position >= l_causal_start[causal_start_counter] and position < l_causal_start[causal_start_counter]+causal_bin_length):
            # need to update v_snp here
            print('[note] next considered variant is also in a causal bin')
            [v_genetic_effect, d_causals] = get_effect_size(position, mac, tau, rho, d_macs_actual, l_macs_actual, v_genetic_effect, v_snp, d_causals, l_flat_s)
        """
        N is number of total individuals
        v_genetic_effect is [Nx1] vector containing effect sizes due to SNPs summed up
        v_snp is [Nx1] vector containing (0,1,2) coding of current SNP under consideration
        z is effect size
        delta is direction of effect
        s is selection coefficient
        rho and tau are inputs from Uricchio model
        l_macs_actual is sorted keys of d_macs_actual
        d_macs_actual is dictionary, where d_macs_actual[mac] = [list of s]; actual means "for this VCF file and this iteration of simulation"
        d_macs is the dictionary of d[mac]=[list of s] from downsampled SFS_CODE simulations
        """

# -----------------------------------------------------
# 7. Determine environmental effect
print('[INFO] Determining environmental effect size')

var_A = np.var(v_genetic_effect)
var_E = var_A * (( 1 - heritability ) / heritability)

v_environment = np.random.normal(0, math.sqrt(var_E), size = np.size(v_genetic_effect)) # sample E from N(0, var_E)
v_phenotype = v_genetic_effect + v_environment

# -----------------------------------------------------
# 8. Determine cutoff for cases.
print('[INFO] Determining case/control cutoff')

v_cutoff = np.copy(v_phenotype)
v_cutoff = np.sort(v_cutoff)
pheno_cutoff = np.percentile(v_cutoff, 100-prevalence)

# -----------------------------------------------------
# 9. Determine l_control_id and l_case_id (list of individual IDs that are controls and cases, respectively)
print('[INFO] Determining l_control_id and l_case_id')

#   DATA:   l_control_id (list, len = n_control)
#           l_case_id (list, len = n_case)

l_control_id = []
l_case_id = []
for i, id_name in enumerate(l_ids):
    if v_phenotype[i] < pheno_cutoff:
        l_control_id.append(id_name)
    else:
        l_case_id.append(id_name)

# -----------------------------------------------------
# 10. Output phenotype data
print('[INFO] Outputting phenotype data')

# causal variants (use d_causals)
causal_variants_filename = pheno_path+'/'+str(base_name)+'.causals.gz'
outfile = gzip.open( causal_variants_filename, 'wb' )
outfile.write('position\teffect_size\n')
for key in d_causals:
    outfile.write(str(key)+'\t'+str(d_causals[key])+'\n')

outfile.close()

# phenotypes for rvtests

# header is fid iid fatid matid sex y1
# format is ID ID 0 0 0 {1,2}, where 1 for control and 2 for case

rvtest_pheno_filename = pheno_path+'/'+str(base_name)+'.rvtests.pheno.gz'
outfile = gzip.open( rvtest_pheno_filename, 'wb' )
outfile.write('fid\tiid\tfatid\tmatid\tsex\ty1\n')
for i, id_name in enumerate(l_control_id):
    outfile.write(id_name+'\t'+id_name+'\t0\t0\t0\t1\n')

for i, id_name in enumerate(l_case_id):
    outfile.write(id_name+'\t'+id_name+'\t0\t0\t0\t2\n')

outfile.close()

# individual genetic, environment, phenotype (gep)
gep_filename = pheno_path+'/'+str(base_name)+'.individual.gep.gz'
outfile = gzip.open( gep_filename, 'wb' )
for i, id_name in enumerate(l_ids):
    outfile.write(id_name+'\t'+str(v_genetic_effect[i])+'\t'+str(v_environment[i])+'\t'+str(v_phenotype[i])+'\n')

outfile.close()

# -----------------------------------------------------
# 11. Generate reference panel ID list and case/control panel ID list
print('[INFO] Generating reference panel ID list and case/control ID list')

# inputs
# l_control_id      list of individual ID's that are controls (all controls in simulated population)
# l_case_id         list of individual ID's that are cases (all cases in simulated population)

# output
# l_ref_id          list of individual ID's that are reference panel individuals (mix of cases and controls, must not take all controls or all cases from simulated population)
# l_cc_id           list of individual ID's that are case/control panel individuals (mix of cases and controls, must roughly match prevalence in population)
# l_remainder_id    list of individual ID's that are in neither reference panel or case/control panel (i.e. held out individuals)

if study_design == 'prop':
    # OPTION 1: random drawing, proportional to population
    # get num_ref individuals from all IDs
    l_cc_id = random.sample(l_ids, num_cc)

    # # get num_cc individuals from all IDs minus the individuals chosen for reference panel
    l_ref_id = random.sample(list(set(l_ids) - set(l_cc_id)), num_ref)

elif study_design == '5050':
    # OPTION 2: half cases and half controls (50/50)
    l_ref_case_id = random.sample(l_case_id, int(math.floor((prevalence/100.)*num_ref)))
    l_ref_control_id = random.sample(l_control_id, int(math.floor(((100-prevalence)/100.)*num_ref)))

    if 0.5*num_cc > len(list(set(l_case_id) - set(l_ref_case_id))):
        print('[ERROR] Insufficient number of cases for 50/50 case/control panel: 0.5*num_cc > len(list(set(l_case_id) - set(l_ref_case_id)))')
    else:
        l_cc_case_id = random.sample(list(set(l_case_id) - set(l_ref_case_id)), int(0.5*num_cc))

    if 0.5*num_cc > len(list(set(l_control_id) - set(l_ref_control_id))):
        print('[ERROR] Insufficient number of controls for 50/50 case/control panel: 0.5*num_cc > len(list(set(l_control_id) - set(l_ref_control_id)))')
    else:
        l_cc_control_id = random.sample(list(set(l_control_id) - set(l_ref_control_id)), int(0.5*num_cc))

    l_ref_id = l_ref_case_id + l_ref_control_id
    l_cc_id = l_cc_case_id + l_cc_control_id

    # ---
    # check reference panel for case/control split
    control = 0
    case = 0
    for id_name in l_ref_id:
        if id_name in l_control_id:
            control = control + 1
        elif id_name in l_case_id:
            case = case + 1
        else:
            print('[ERROR]: ID is not a case or control!')

    print('controls in ref panel: '+str(control))
    print('cases in ref panel: '+str(case))
    print('proportion of cases in ref panel: '+str(float(case)/(case+control)))

    # check case/control panel for case/control split
    control = 0
    case = 0
    for id_name in l_cc_id:
        if id_name in l_control_id:
            control = control + 1
        elif id_name in l_case_id:
            case = case + 1
        else:
            print('[ERROR]: ID is not a case or control!')

    print('controls in cc panel: '+str(control))
    print('cases in cc panel: '+str(case))
    print('proportion of cases in cc panel: '+str(float(case)/(case+control)))


elif study_design == 'extremes':
    # OPTION 3: Extremes
    d_pheno = defaultdict(float)
    for i, id_name in enumerate(l_ids):
        d_pheno[id_name]=v_phenotype[i]

    len([k for k, v in d_pheno.items() if v >= pheno_cutoff])

    tes = d_pheno.values()
    tes.sort()

    extreme_case_cutoff = np.percentile(tes, 95)
    extreme_control_cutoff = np.percentile(tes, 5)

    l_extreme_cases = [k for k, v in d_pheno.items() if v >= extreme_case_cutoff]
    l_extreme_controls = [k for k, v in d_pheno.items() if v <= extreme_control_cutoff]

    print('number of extreme cases: '+str(len(l_extreme_cases)))
    print('number of extreme controls: '+str(len(l_extreme_controls)))

    l_cc_id = l_extreme_cases+l_extreme_controls
    l_ref_id = random.sample(list(set(l_ids) - set(l_cc_id)), num_ref)
else:
    print('[ERROR] Study design not recognized')

# -----------------------------------------------------
# 12. Split input VCF into reference panel VCF
# case/control panel and remainder panel VCF's come after downsampling
print('[INFO] Extracting reference panel VCF from input VCF')

# extract reference panel individuals to ref.vcf.gz
ref_sample_file_name = rvtest_input_path+'/'+base_name+'.ref.sample_file'
ref_vcf_gz = rvtest_input_path+'/'+base_name+'.ref.vcf.gz'

outfile = open(ref_sample_file_name, 'w')

for id_name in l_ref_id:
    outfile.write(id_name+'\n')

outfile.close()

bcf_view_call = '/netapp/home/dmctong/programs/bcftools/bcftools view '
dash_sample_file = '--samples-file '+str(ref_sample_file_name)+' '
dash_o = '-o ' + str(ref_vcf_gz) + ' --output-type z '
dash_i = vcf_filename
chunk_call = bcf_view_call + dash_sample_file + dash_o + dash_i
subprocess.call(chunk_call, shell=True)

def tabix(input_vcf_gz):
    tabix_call='/netapp/home/dmctong/programs/htslib/tabix -fp vcf '+input_vcf_gz
    subprocess.call(tabix_call,shell=True)

tabix(ref_vcf_gz)

# -----------------------------------------------------
# 13. Downsample cc.vcf to cc.ds.vcf (by matching genotyping array)
print('[INFO] Downsampling case/control VCF to match genotyping array')

# read in SFS for a particular genotyping array; currently working with Illumina Omni2.5-8
# RETURNS: d_maf_histogram[(min_MAF_in_bin,max_MAF_in_bin)]=count_in_bin
def read_in_sfs_from_genotyping_array(array_sfs_file):
    array_sfs=open(array_sfs_file,'r')
    d_maf_histogram=defaultdict(int)
    for line in array_sfs:
        if line.startswith('MAF'):
            continue
        if line.startswith('0'):
            continue
        line=line.strip().split()
        d_maf_histogram[(float(line[0][1:-1]),float(line[1][2:]))]=int(line[5])
    return d_maf_histogram

# read in SFS from genotyping array data file (Illumina Omni2.5 in this case)
d_maf_histogram = read_in_sfs_from_genotyping_array(array_sfs_file)

target_distance = 1230 # INPUT

d_garray = defaultdict(list)
for i,key in enumerate(d_maf_histogram):
    d_garray[i] = d_maf_histogram[key]

maf_cutoffs = np.asarray([0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5])
mac_cutoffs = maf_cutoffs*(2*ns)
d_vcf_histogram = defaultdict(list)
for mac in d_sfs:
    for i,cutoff in enumerate(mac_cutoffs):
        if mac<cutoff:
            #print(mac,cutoff)
            d_vcf_histogram[i-1] = d_vcf_histogram[i-1]+d_sfs[mac]
            break

# for key in d_vcf_histogram:
#     print(key,len(d_vcf_histogram[key]))
l_positions = []
l_prob = np.asarray(d_garray.values())/float(sum(d_garray.values()))
j = 0
# region_length = max(d_sfs.values())
mean_distance = region_length
while mean_distance > target_distance or j < 1000:
    maf_bin = np.random.choice(d_garray.keys(), size=1, p=l_prob)
    #random_position = np.random.choice(d_vcf_histogram[maf_bin[0]],size = 1)
    l_test = d_vcf_histogram[maf_bin[0]]
    if not l_test:
        continue
    random_position = l_test[np.random.randint(0, len(l_test))]
    d_vcf_histogram[maf_bin[0]].remove(random_position)
    l_positions.append(random_position)
    l_positions.sort()
    l_distance = [0]*len(l_positions)
    for i,thing in enumerate(l_positions):
        if i==0 or i==len(l_positions):
            continue
        l_distance[i] = l_positions[i] - l_positions[i-1]
    j = j+1
    mean_distance = np.mean(l_distance[1:-1])
    if j%100==0:
        print(j, mean_distance)

l_positions.sort()

# Debugging
# print('[INFO] Printing out l_positions, l_distance, mean_distance, target_distance, j')
# print('l_positions')
# print(l_positions)
# print('l_distance')
# print(l_distance)
# print('mean_distance')
# print(mean_distance)
# print('target_distance')
# print(target_distance)
# print('j')
# print(j)


# -----------------------------------------------------
# 14. Create case/control panel VCF
print('[INFO] Creating case/control panel VCF')

cc_vcf_gz = rvtest_input_path+'/'+base_name+'.cc.vcf.gz'
cc_ds_vcf_gz = rvtest_input_path+'/'+base_name+'.cc.ds.vcf.gz'

# sample file - defines which samples are in cc panel
cc_sample_filename = rvtest_input_path+'/'+base_name+'.cc.sample_file'
outfile = open(cc_sample_filename, 'w')
for id_name in l_cc_id:
    outfile.write(id_name+'\n')

outfile.close()

# regions file - defines positions in genotyping array
regions_filename = rvtest_input_path+'/'+base_name+'.genotyped_positions.txt'
outfile=open(regions_filename,'w')
for thing in l_positions:
    outfile.write('1\t'+str(thing)+'\n')

outfile.close()

# cc_vcf_gz for full data
bcf_view_call = '/netapp/home/dmctong/programs/bcftools/bcftools view '
dash_sample_file = '--samples-file '+str(cc_sample_filename)+' '
dash_o = '-o ' + str(cc_vcf_gz) + ' --output-type z '
dash_i = bgz_vcf_filename # needs to be bgzip and tabix??!
chunk_call = bcf_view_call + dash_sample_file + dash_o + dash_i
subprocess.call(chunk_call, shell=True)

tabix(cc_vcf_gz)

# reduce input VCF to cc panel individuals and genotyping array positions (cc_ds_vcf_gz)
bcf_view_call = '/netapp/home/dmctong/programs/bcftools/bcftools view '
dash_regions_file = '--regions-file '+str(regions_filename)+' ' # NB regions_file from downsampling to genotyping array
dash_o = '-o ' + str(cc_ds_vcf_gz) + ' --output-type z '
dash_i = cc_vcf_gz # needs to be bgzip and tabix??!
chunk_call = bcf_view_call + dash_regions_file + dash_o + dash_i
subprocess.call(chunk_call, shell=True)

tabix(cc_ds_vcf_gz)

# -----------------------------------------------------
# 15. Conversions
# convert cc.ds.vcf.gz to cc.ds.gen and cc.ds.gen.sample for SHAPEIT/IMPUTE4
# convert ref.vcf.gz into haps/legend/sample format

print('[INFO] Converting case/control VCF to .gen files')

def convert_vcf_to_gen(input_vcf, gen_file, gen_sample_file): # gen and gen_sample files are output, VCF file is input
    conversion_call = '/netapp/home/dmctong/programs/bcftools/bcftools convert --gensample '+gen_file+','+gen_sample_file+' '+input_vcf
    subprocess.call(conversion_call, shell=True)

cc_ds_gen = rvtest_input_path+'/'+base_name+'.cc.ds.gen'
cc_ds_gen_sample = rvtest_input_path+'/'+base_name+'.cc.ds.gen.sample'
convert_vcf_to_gen(cc_ds_vcf_gz, cc_ds_gen, cc_ds_gen_sample)

print('[INFO] Converting ref VCF to haplegendsample files')

def convert_vcf_to_hls(input_vcf, haps_file, legend_file, sample_file): # HLS files are output!
    conversion_call='/netapp/home/dmctong/programs/bcftools/bcftools convert --haplegendsample '+haps_file+','+legend_file+','+sample_file+' '+input_vcf
    subprocess.call(conversion_call, shell=True)

ref_haps = rvtest_input_path+'/'+base_name+'.ref.haps'
ref_legend = rvtest_input_path+'/'+base_name+'.ref.legend'
ref_sample = rvtest_input_path+'/'+base_name+'.ref.sample'
convert_vcf_to_hls(ref_vcf_gz, ref_haps, ref_legend, ref_sample)

subprocess.call('gzip '+ref_haps, shell=True)
subprocess.call('gzip '+ref_legend, shell=True)
subprocess.call('gzip '+ref_sample, shell=True)

ref_haps_gz = ref_haps+'.gz'
ref_sample_gz = ref_sample+'.gz'
ref_legend_gz = ref_legend+'.gz'

# Testing the VCF itself
# bcftools stats [OPTIONS] A.vcf.gz [B.vcf.gz]
bcf_stats_call = '/netapp/home/dmctong/programs/bcftools/bcftools stats -v '
bash_call = bcf_stats_call+bgz_vcf_filename
subprocess.call(bash_call, shell=True)

# SFS

# Delete ref_vcf_gz since it is unneeded in P3
subprocess.call('rm '+ref_vcf_gz,shell=True)