# run this file in iqint
# it is more efficient to do double for loops than wait in line?

import subprocess
from collections import defaultdict
import os
import numpy as np
import sys

# from subprocess import Popen, PIPE
# unique_id = 'EUR-10k-5050'
# vcf_file = 'EUR-20000.5mb.vcf.gz'
# num_sims = 100
# num_params = 16

# ${VCF_FILENAME_SHORT} ${UNIQUE_ID} ${NUM_SIMS} ${NUM_PARAMS}
vcf_file = str(sys.argv[1])
unique_id = str(sys.argv[2])
num_sims = int(sys.argv[3])
num_params = int(sys.argv[4])

d_results=defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(list)))))
d_allsims = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(float))))) # d_allsims[param_num][ve][num_correct]; d_allsims[param_num][ve][num_total]

mkdir_call='mkdir /hernandez/mandrill/users/dmctong/imp_rvat/data/'+str(unique_id)+'/cat_data'
print(mkdir_call)
subprocess.call(mkdir_call,shell=True)

for param_num in range(1,num_params+1):
    for iter_num in range(1,num_sims+1): 
        base_name = str(param_num)+'.'+str(iter_num)
        print(base_name)
        # directories
        rvtest_output_dir = '/hernandez/mandrill/users/dmctong/imp_rvat/data/'+str(unique_id)+'/'+str(param_num)+'/'+str(iter_num)+'/rvtest_output'
        gwas_output_dir = '/hernandez/mandrill/users/dmctong/imp_rvat/data/'+str(unique_id)+'/'+str(param_num)+'/'+str(iter_num)+'/gwas'
        # cat rvtest outputs
        rvtest_cat_call = 'cat '+rvtest_output_dir+'/*.rvtest.data.all > '+rvtest_output_dir+'/'+str(base_name)+'.cat.all'
        subprocess.call(rvtest_cat_call, shell=True)
        # find unique lines only (repeats on full data)
        awk_call = 'awk \'!seen[$0]++\' '+rvtest_output_dir+'/'+base_name+'.cat.all'+' > '+rvtest_output_dir+'/'+base_name+'.cat.all.uniq'
        subprocess.call(awk_call, shell=True)
        # remove non-unique rvtest results
        rm_call = 'rm '+rvtest_output_dir+'/'+base_name+'.cat.all'
        subprocess.call(rm_call, shell=True)
        # cat firth results
        firth_cat_call = 'cat '+rvtest_output_dir+'/*.firth.data.all >> '+rvtest_output_dir+'/'+str(base_name)+'.cat.all.uniq'
        subprocess.call(firth_cat_call, shell=True)
        # cat gwas outputs
        gwas_cat_call = 'cat '+gwas_output_dir+'/*.gwas.data.all >> '+rvtest_output_dir+'/'+str(base_name)+'.cat.all.uniq'
        subprocess.call(gwas_cat_call, shell=True)
        # populate d_results
        file_name = '/hernandez/mandrill/users/dmctong/imp_rvat/data/'+str(unique_id)+'/'+str(param_num)+'/'+str(iter_num)+'/rvtest_output/'+str(base_name)+'.cat.all.uniq'
        if os.path.isfile(file_name):
            infile = open(file_name,'r')
            for x,line in enumerate(infile):
                line = line.strip().split()
                test = str(line[0])
                imp = int(line[1])
                ve = float(line[7])
                rvtest_result = int(float(line[8]))
                d_results[iter_num][test][imp][ve][param_num].append(rvtest_result)
                # pseudocode for a moment here:
                # if ve!=0 and rvtest_result == 1: it is correct
                # if ve!=0 and rvtest_result == 0: it is wrong
                # if ve==0 and rvtest_result == 1: it is wrong
                # if ve==0 and rvtest_result == 0: it is correct
                d_allsims[param_num][test][imp][ve]['num_total'] = d_allsims[param_num][test][imp][ve]['num_total'] + 1
                if ve != 0:
                    if rvtest_result == 1:
                        d_allsims[param_num][test][imp][ve]['num_correct'] = d_allsims[param_num][test][imp][ve]['num_correct'] + 1
                elif ve == 0:
                    if rvtest_result == 0:
                        d_allsims[param_num][test][imp][ve]['num_correct'] = d_allsims[param_num][test][imp][ve]['num_correct'] + 1
        else:
            print('[ERROR] File does not exist: '+str(file_name))


# STEP 2: Calculate mean point
# d_allsims = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(float))))) # d_allsims[param_num][ve][num_correct]; d_allsims[param_num][ve][num_total]
# outfile=open('/hernandez/mandrill/users/dmctong/imp_rvat/data/'+str(unique_id)+'/cat_data/results_bootstrap.txt','w')
# for iter_num in d_results:
#     for test in d_results[iter_num]:
#         for imp in d_results[iter_num][test]:
#             for ve in d_results[iter_num][test][imp]:
#                 for param_num in d_results[iter_num][test][imp][ve]:
#                     if len(d_results[iter_num][test][imp][ve][param_num])>0:
#                         if ve==0:
#                             numerator = float(len(d_results[iter_num][test][imp][ve][param_num]))-float(sum(d_results[iter_num][test][imp][ve][param_num]))
#                             d_allsims[param_num][test][imp][ve]['num_correct'] = d_allsims[param_num][test][imp][ve]['num_correct'] + numerator
#                             d_allsims[param_num][test][imp][ve]['num_total']=d_allsims[param_num][test][imp][ve]['num_total']+float(len(d_results[iter_num][test][imp][ve][param_num]))
#                         else:
#                             d_allsims[param_num][test][imp][ve]['num_correct']=d_allsims[param_num][test][imp][ve]['num_correct']+float(sum(d_results[iter_num][test][imp][ve][param_num]))
#                             d_allsims[param_num][test][imp][ve]['num_total']=d_allsims[param_num][test][imp][ve]['num_total']+float(len(d_results[iter_num][test][imp][ve][param_num]))

# # for quick checking:
# for param_num in d_allsims:
#     for test in d_allsims[param_num]:
#         for imp in d_allsims[param_num][test]:
#             for ve in d_allsims[param_num][test][imp]:
#                 # for num in d_allsims[param_num][test][imp][ve]:
#                 print(param_num,test,imp,ve,d_allsims[param_num][test][imp][ve])

# STEP 3: bootstrap to find error bars

# d_results[iter_num][test][imp][ve][param_num]=[rvtest_result]

# sample 100 sims with replacement 1000 times
# repeat 1000 times:
# pick 100 random sims with replacement
# calculate num_correct/num_total

# d_bootstrap[param_num][test][imp][ve]=[bootstrapped power]
d_bootstrap=defaultdict(lambda: defaultdict(lambda:defaultdict(lambda: defaultdict(list))))
for i in range(0,100):
    print('bootstrap iteration: '+str(i))
    l_random_sims=np.random.choice(range(1,num_sims+1),size=100)
    for sim in l_random_sims:
        for test in d_results[sim]:
            for imp in d_results[sim][test]:
                for ve in d_results[sim][test][imp]:
                    for param_num in d_results[sim][test][imp][ve]:
                        if len(d_results[sim][test][imp][ve][param_num])>0:
                            if ve==0:
                                correct = len(d_results[sim][test][imp][ve][param_num])-sum(d_results[sim][test][imp][ve][param_num])
                            else:
                                correct = sum(d_results[sim][test][imp][ve][param_num])
                            total = len(d_results[sim][test][imp][ve][param_num])
                            power = float(correct)/float(total)
                            d_bootstrap[param_num][test][imp][ve].append(power)

d_bootstrap_results=defaultdict(lambda: defaultdict(lambda:defaultdict(lambda: defaultdict(lambda: defaultdict(float)))))
for test in d_results[sim]:
    for imp in d_results[sim][test]:
        for ve in d_results[sim][test][imp]:
            for param_num in d_results[sim][test][imp][ve]:
                if len(d_results[sim][test][imp][ve][param_num])>0:
                    lower = np.percentile(d_bootstrap[param_num][test][imp][ve],25)
                    upper = np.percentile(d_bootstrap[param_num][test][imp][ve],75)
                    d_bootstrap_results[param_num][test][imp][ve]['low'] = lower
                    d_bootstrap_results[param_num][test][imp][ve]['high'] = upper

# STEP 4: OUTPUT
outfile=open('/hernandez/mandrill/users/dmctong/imp_rvat/data/'+str(unique_id)+'/cat_data/results_bootstrap.txt','w')
for param_num in d_allsims:
    for test in d_allsims[param_num]:
        for imp in d_allsims[param_num][test]:
            for ve in d_allsims[param_num][test][imp]:
                # for num in d_allsims[param_num][test][imp][ve]:
                print(param_num,test,imp,ve,d_allsims[param_num][test][imp][ve])
                num_correct = d_allsims[param_num][test][imp][ve]['num_correct']
                num_total = d_allsims[param_num][test][imp][ve]['num_total']
                mean = float(num_correct)/float(num_total)
                low = d_bootstrap_results[param_num][test][imp][ve]['low']
                high = d_bootstrap_results[param_num][test][imp][ve]['high']
                outfile.write(str(param_num)+'\t'+str(test)+'\t'+str(imp)+'\t'+str(ve)+'\t'+str(num_correct)+'\t'+str(num_total)+'\t'+str(mean)+'\t'+str(low)+'\t'+str(high)+'\n')

outfile.close()

# rsync /hernandez/mandrill/users/dmctong/imp_rvat/data/2018-04-30-1/cat_data/results_bootstrap.txt dominic@169.230.84.121:/home/dominic/Dropbox/lab/imp_rvat/2018-04-30-1