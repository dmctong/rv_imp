# combine get_5mb_from_vcf and initial phenotype simulation file-making

# Extract 5Mb from simulated VCF
# gzip while keeping original
# bgzip into ./input-bgz (this will overwrite original .gz.....)
# tabix into ./input-bgz
# copy gzip to scrapp
# copy bgzip+tabix to scrapp

# get all input variables
# get paths
# generate ad_filename and vcf_mac_filename

# NB YOU STILL NEED TO RSYNC VCF.GZ AND TABIX FILES TO SCRAPP FOR USE
# rsync /hernandez/mandrill/users/dmctong/msprime/data/input-gzip/EUR-20000.5mb.vcf.gz /scrapp/dmctong/msprime/data/input-gzip
# rsync /hernandez/mandrill/users/dmctong/msprime/data/input-bgz/EUR-20000.5mb.vcf.gz* /scrapp/dmctong/msprime/data/input-bgz

import gzip
import subprocess
import os.path
import sys

# ===============================================
# Extract 5Mb region from chr22 simulated file
# ===============================================

# the msprime output is expected in input_gz_path, so rsync it there!

#-------------------------------------
# INPUTS
p1_jobid = str(sys.argv[1])
population = str(sys.argv[2])
pop_size = str(sys.argv[3])

start_pos = 17000000
end_pos = 22000000

vcf_file_gz = str(p1_jobid)+'.'+str(population)+'-'+str(pop_size)+'.vcf.gz'
vcf_filename_short = str(p1_jobid)+'.'+str(population)+'-'+str(pop_size)+'.5mb.vcf.gz'

# vcf_file_gz = 'EUR-10000.vcf.gz'
# vcf_filename_short = 'EUR-10000.5mb.vcf.gz' # used in next section (generate mac and ad files)

input_gz_path = '/hernandez/mandrill/users/dmctong/imp_rvat/data/input-gzip' # P1 generates output to data

msprime_output_vcf_gz = input_gz_path+'/'+vcf_file_gz
vcf_5mb_unzip = msprime_output_vcf_gz[:-7]+'.5mb.vcf'
vcf_5mb_gz = msprime_output_vcf_gz[:-7]+'.5mb.vcf.gz'

input_bgz_path = '/hernandez/mandrill/users/dmctong/imp_rvat/data/input-bgz'
vcf_5mb_pre_bgz = input_bgz_path+'/'+vcf_file_gz[:-7]+'.5mb.vcf'
vcf_5mb_bgz = input_bgz_path+'/'+vcf_file_gz[:-7]+'.5mb.vcf.gz'

infile = gzip.open(msprime_output_vcf_gz,'rb')
outfile = open(vcf_5mb_unzip, 'w')

for i,line in enumerate(infile):
    if line.startswith('##'):
        outfile.write(line)
        continue
    elif line.startswith('#CHROM'):
        outfile.write(line)
        continue
    else:
        line2=line.strip().split()
        pos=int(line2[1])
        if pos < start_pos:
            continue
        elif start_pos <= pos <= end_pos:
            outfile.write(line)
        elif pos>end_pos:
            print('[INFO] pos > end_pos')
            break
    if i%100000==0:
        print(i, pos)

infile.close()
outfile.close()

# next steps: currently have .vcf
# 1. gzip while keeping original
# 2. bgzip into ./input-bgz (this will overwrite original .gz.....)
# 3. tabix into ./input-bgz
# 4. copy gzip to scrapp
# 5. copy bgzip+tabix to scrapp

# gzip -c /hernandez/mandrill/users/dmctong/msprime/data/EUR-20000.5mb.vcf > /hernandez/mandrill/users/dmctong/msprime/data/input-gzip/EUR-20000.5mb.vcf.gz

subprocess.call('gzip -c '+vcf_5mb_unzip+' > '+vcf_5mb_gz, shell=True)

# mv vcf_5mb_unzip to ./input-bgz/
subprocess.call('mv '+vcf_5mb_unzip+' '+input_bgz_path, shell=True)

subprocess.call('/netapp/home/dmctong/programs/htslib/bgzip -f '+vcf_5mb_pre_bgz,shell=True)
subprocess.call('/netapp/home/dmctong/programs/htslib/tabix -fp vcf '+vcf_5mb_bgz,shell=True)

# ===============================================
# Generate vcf_mac_filename and ad_filename for phenotype simulation
# ===============================================

#-------------------------------------
# Paths and filenames

mandrill_path = '/hernandez/mandrill/users/dmctong/imp_rvat/data'
support_path = mandrill_path+'/support'

vcf_filename = mandrill_path+'/input-gzip/'+str(vcf_filename_short) #TODO rsync from here to temp_dir?

ad_filename = support_path+'/'+vcf_filename_short[:-7]+'.ad.even.gz' #additive model file name
vcf_mac_filename = support_path+'/'+vcf_filename_short[:-7]+'.mac.gz'

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
                    if i <= 8:
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

ad_filename_odd = support_path+'/'+vcf_filename_short[:-7]+'.ad.odd.gz' #additive model file name
subprocess.call('cp '+str(ad_filename)+' '+str(ad_filename_odd),shell=True) # create even and odd version to reduce read/write conflicts