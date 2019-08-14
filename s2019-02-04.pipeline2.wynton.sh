#$ -S /bin/sh
#$ -cwd
#$ -r y
#$ -o /scrapp/dmctong/imp_rvat/50k_scripts/2019-04/
#$ -e /scrapp/dmctong/imp_rvat/50k_scripts/2019-04/
#$ -l mem_free=5G
#$ -l scratch=5G
#$ -l h_rt=336:00:00
#$ -t 2-64:2

# ================================================
# INPUT VARIABLES
# ================================================
# VCF_FILENAME_SHORT=EUR-20000.5mb.vcf.gz

ITERATION_NUM=$1 # 1-100, denotes replicates
P1_JOBID=$2 # job ID of P1 run; VCF file is named ${P1_JOBID}.${POPULATION}-${POP_SIZE}
POPULATION=$3
POP_SIZE=$4
STUDY_DESIGN=$5

UNIQUE_ID=${POPULATION}-${POP_SIZE}-${STUDY_DESIGN}-3

EXIT_STATUS=0

if [[ $EXIT_STATUS -eq "0" ]]; then
        TEMP_DIR=`mktemp -d -p /scratch/`
        EXIT_STATUS=$(( ${EXIT_STATUS} || $? ))
        if [[ $EXIT_STATUS -ne "0" ]]; then
                echo "FAIL: Making working directory failed"
        fi
fi

# ================================================
# JOB VARIABLES
# ================================================

START=$(date +%s)
echo start_time: $START
echo SGE_TASK_ID: $SGE_TASK_ID
echo JOB_ID: $JOB_ID
echo TEMP_DIR: $TEMP_DIR

if [[ $EXIT_STATUS -eq "0" ]]; then
        read REGION_START REGION_END NS NB CAUSAL_BIN_LENGTH NUM_CAUSAL_BINS RHO TAU HERITABILITY PREVALENCE NUM_REF NUM_CC <<< $(sed "${SGE_TASK_ID}q;d" /netapp/home/dmctong/imp_rvat/50k_scripts/${UNIQUE_ID}.pipeline2.todo.txt)
        EXIT_STATUS=$(( ${EXIT_STATUS} || $? ))
        if [[ $EXIT_STATUS -ne "0" ]]; then
                echo "FAIL: reading todo failed"
        fi
fi

echo REGION_START: ${REGION_START}
echo REGION_END: ${REGION_END}
echo NS: ${NS}
echo NB: ${NB}
echo CAUSAL_BIN_LENGTH: ${CAUSAL_BIN_LENGTH}
echo NUM_CAUSAL_BINS: ${NUM_CAUSAL_BINS}
echo RHO: ${RHO}
echo TAU: ${TAU}
echo HERITABILITY: ${HERITABILITY}
echo PREVALENCE: ${PREVALENCE}
echo NUM_REF: ${NUM_REF}
echo NUM_CC: ${NUM_CC}
echo SGE_TASK_ID: ${SGE_TASK_ID}
echo TEMP_DIR: ${TEMP_DIR}
echo ITERATION_NUM: ${ITERATION_NUM}
echo P1_JOBID: ${P1_JOBID}
echo POPULATION: ${POPULATION}
echo POP_SIZE: ${POP_SIZE}
echo STUDY_DESIGN: ${STUDY_DESIGN}


# echo START_INDEX: ${START_INDEX} #P3 and P4 only
# echo END_INDEX: ${END_INDEX} #P3 and P4 only
# echo RVTEST_SET_LENGTH: ${RVTEST_SET_LENGTH} #P3 and P4 only

# scl enable python27 bash
# LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/$HOME/gsl-2.4/lib
# L_POP=( 1000 2000 5000 10000 20000 50000 100000 200000 )
# POP_SIZE=${L_POP[i]}

# ---------------------
# job start rate limiter
# ---------------------
JOB_START_IO=/netapp/home/dmctong/job_start # presence of a file in this dir denotes actively running job
N_JOB_START=20 # Number of parallel jobs running allowed
# Sleep until number of jobs started falls below specified level
while [ `ls -1 $JOB_START_IO | wc -l ` -ge $N_JOB_START ]; do sleep 15; done
# record job start
touch ${JOB_START_IO}/${USER}_${JOB_ID}_${SGE_TASK_ID}
rm ${JOB_START_IO}/${USER}_${JOB_ID}_${SGE_TASK_ID}

# ---------------------
# Do the job
# ---------------------
if [[ $EXIT_STATUS -eq "0" ]]; then
        python /netapp/home/dmctong/imp_rvat/50k_scripts/s2019-02-04.pipeline2.wynton.py ${REGION_START} ${REGION_END} ${NS} ${NB} ${CAUSAL_BIN_LENGTH} ${NUM_CAUSAL_BINS} ${RHO} ${TAU} ${HERITABILITY} ${PREVALENCE} ${NUM_REF} ${NUM_CC} ${SGE_TASK_ID} ${TEMP_DIR} ${ITERATION_NUM} ${P1_JOBID} ${POPULATION} ${POP_SIZE} ${STUDY_DESIGN}
        # TODO change inputs and script name
        EXIT_STATUS=$(( ${EXIT_STATUS} || $? ))
        if [[ $EXIT_STATUS -ne "0" ]]; then
                echo "FAIL: pipeline2.py failed"
        fi
fi

# ---------------------
# read write rate limiter
# ---------------------
# mkdir /hernandez/mandrill/users/dmctong/imp_rvat/data/${UNIQUE_ID}
# mkdir /hernandez/mandrill/users/dmctong/imp_rvat/data/${UNIQUE_ID}/${SGE_TASK_ID}
# mkdir /hernandez/mandrill/users/dmctong/imp_rvat/data/${UNIQUE_ID}/${SGE_TASK_ID}/${ITERATION_NUM}

mkdir /wynton/scratch/dmctong/imp_rvat/data/${UNIQUE_ID}
mkdir /wynton/scratch/dmctong/imp_rvat/data/${UNIQUE_ID}/${SGE_TASK_ID}
mkdir /wynton/scratch/dmctong/imp_rvat/data/${UNIQUE_ID}/${SGE_TASK_ID}/${ITERATION_NUM}

MANDRILL_IO=/netapp/home/dmctong/mandrill_io # The presence of a file in this directory denotes an active read/write
N_READ_WRITE=20 # Number of parallel read/writes allowed
#Sleep until number of read/writes to mandrill drops below specified level
while [ `ls -1 $MANDRILL_IO | wc -l ` -ge $N_READ_WRITE ]; do sleep 15; done
#Record IO start
touch ${MANDRILL_IO}/${USER}_${JOB_ID}_${SGE_TASK_ID}
#Perform I/O to/from mandrill

if [[ $EXIT_STATUS -eq "0" ]]; then
        rsync -rauz ${TEMP_DIR}/* /wynton/scratch/dmctong/imp_rvat/data/${UNIQUE_ID}/${SGE_TASK_ID}/${ITERATION_NUM}
        EXIT_STATUS=$(( ${EXIT_STATUS} || $? ))
        if [[ $EXIT_STATUS -ne "0" ]]; then
                echo "FAIL: Return transfer failed"
        fi
fi

#Record IO end
rm ${MANDRILL_IO}/${USER}_${JOB_ID}_${SGE_TASK_ID}

rm -r ${TEMP_DIR}

qstat -j ${JOB_ID}

END=$(date +%s)
echo end_time: $END
echo Run time was "$(expr $END - $START)" seconds

if [ $EXIT_STATUS -eq "0" ]; then
        echo "Complete"
else
        echo "Complete with Errors"
fi