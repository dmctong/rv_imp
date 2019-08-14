import msprime
import gzip
import os
import math
import time
import sys

#-------------------------------------
# Inputs
diploid_eur_sample_size = int(sys.argv[1])
job_id = str(sys.argv[2])
outfilename = str(sys.argv[3]) # should be to $UNIQUE_NAME.vcf.gz
temp_dir = str(sys.argv[4]) # so outfile goes to temp_dir/outfilename

outfilename_with_dir = temp_dir+'/'+outfilename

eur_sample_size = diploid_eur_sample_size*2 # number of haploid alleles to be simulated; true sample size is eur_sample_size/2

#-------------------------------------
# Read in recombination rate file
rr_file = "/wynton/scratch/dmctong/imp_rvat/data/genetic_map_chr22_combined_b37.20140701.add_chr22.txt"
recomb_map = msprime.RecombinationMap.read_hapmap(rr_file)

#-------------------------------------
# Define Tennessen model in msprime
def out_of_africa(eur_sample_size):
    # First we set out the maximum likelihood values of the various parameters
    # given in Table 1.
    N_A = 7300
    N_B = 1861
    N_AF0 = 14474
    N_EU0 = 1032
    N_AS0 = 550
    # Times are provided in years, so we convert into generations.
    generation_time = 25
    T_AF = 148e3 / generation_time
    T_B = 51e3 / generation_time
    T_EU_AS = 23e3 / generation_time
    T_EXP = 5.1e3 / generation_time
    # We need to work out the starting (diploid) population sizes based on
    # the growth rates provided for these two populations
    r_EU_not_exp = 0.00307
    r_AS = 0.0048
    r_EU_exp = 0.0195
    r_AFR_exp = 0.0166
    N_EU = (N_EU0 / math.exp(-r_EU_not_exp * (T_EU_AS-T_EXP))) / math.exp(-r_EU_exp * T_EXP)
    N_AS = (N_AS0 / math.exp(-r_AS * T_EU_AS))
    N_AF = (N_AF0 / math.exp(-r_AFR_exp * (T_EXP)))
    # Migration rates during the various epochs.
    m_AF_B = 15e-5
    m_AF_EU = 2.5e-5
    m_AF_AS = 0.78e-5
    m_EU_AS = 3.11e-5
    # Population IDs correspond to their indexes in the population
    # configuration array. Therefore, we have 0=YRI, 1=CEU and 2=CHB
    # initially.
    population_configurations = [
        # African population has size N_AFR and growth rate r_AFR_exponential
        msprime.PopulationConfiguration(
            sample_size=0, initial_size=N_AF, growth_rate=r_AFR_exp),
        # European population has size N_EUR and growth rate r_EUR_exponential
        msprime.PopulationConfiguration(
            sample_size=eur_sample_size, initial_size=N_EU, growth_rate=r_EU_exp),
        # Asian population has size N_ASN and growth rate r_ASN
        msprime.PopulationConfiguration(
            sample_size=0, initial_size=N_AS, growth_rate=r_AS)
    ]
    migration_matrix = [
        [      0, m_AF_EU, m_AF_AS],
        [m_AF_EU,       0, m_EU_AS],
        [m_AF_AS, m_EU_AS,       0],
    ]
    demographic_events = [
        #-------------------------------------
        # At 5.1kya (T_exp):
        # European population changes growth rate to r_EUR_not_exponential
        # NB if initial_size=None, population size is calculated from previous growth rate/size
        msprime.PopulationParametersChange(
            time=T_EXP, initial_size=None, growth_rate=r_EU_not_exp, population_id=1),
        # African population changes growth rate to 0
        msprime.PopulationParametersChange(
            time=T_EXP, initial_size=None, growth_rate=0, population_id=0),
        #-------------------------------------
        # At 23 kya (T_EU_AS):
        # EUR has size 1032
        # ASN has size 550
        # EUR and ASN merge to form B
        msprime.MassMigration(
            time=T_EU_AS, source=2, destination=1, proportion=1.0),
        # B has initial size N_B and no growth rate (check this)
        msprime.PopulationParametersChange(
            time=T_EU_AS, initial_size=N_B, growth_rate=0, population_id=1),
        # B has migration rate with AFR
        msprime.MigrationRateChange(
            time=T_EU_AS, rate=m_AF_B, matrix_index=(1, 0)),
        # AFR has migration rate with B
        msprime.MigrationRateChange(
            time=T_EU_AS, rate=m_AF_B, matrix_index=(0, 1)),
        #-------------------------------------
        # At 51kya (T_B):
        # B has size 1861
        # B merges into AFR
        msprime.MassMigration(
            time=T_B, source=1, destination=0, proportion=1.0),
        # AFR has size 14474
        #-------------------------------------
        # At 148kya (T_AF):
        # AFR has size 7300
        msprime.PopulationParametersChange(
            time=T_AF, initial_size=N_A, population_id=0)
    ]
    # Use the demography debugger to print out the demographic history
    # that we have just described.
    dp = msprime.DemographyDebugger(
        Ne=N_A,
        population_configurations=population_configurations,
        migration_matrix=migration_matrix,
        demographic_events=demographic_events)
    dp.print_history()
    tree_sequence = msprime.simulate(recombination_map=recomb_map,
        mutation_rate=1e-8,
        population_configurations=population_configurations,
        migration_matrix=migration_matrix,
        demographic_events=demographic_events)
    return tree_sequence

#-------------------------------------
# simulate the sequences
start_time=time.time()
tree_sequence=out_of_africa(eur_sample_size)
end_time=time.time()
sim_time = end_time - start_time
print('simulation duration (s): '+str(sim_time))

#-------------------------------------
# output to VCF
start_time=time.time()
vcf_file=gzip.open(outfilename_with_dir,'wb')
try:
    tree_sequence.write_vcf(vcf_file,2)
finally:
    vcf_file.close()

end_time=time.time()
write_time = end_time - start_time
print('writing to outfile time (s): '+str(write_time))
print('Total runtime (s): '+str(write_time + sim_time))