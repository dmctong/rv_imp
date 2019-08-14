#$ -S /bin/sh
#$ -cwd
#$ -r y
#$ -o /scrapp/dmctong/imp_rvat/50k_scripts/2019-03/
#$ -e /scrapp/dmctong/imp_rvat/50k_scripts/2019-03/
#$ -l netapp=2G
#$ -l mem_free=2G
#$ -l scratch=5G
#$ -l h_rt=336:00:00
#$ -l arch=linux-x64
#$ -t 1

# ================================================
# INPUT VARIABLES
# ================================================
START=$(date +%s)
echo start_time: $START

NUM_SIMS=$1
NUM_PARAMS=$2
P1_JOBID=$3
POPULATION=$4
POP_SIZE=$5
STUDY_DESIGN=$6

UNIQUE_ID=${POPULATION}-${POP_SIZE}-${STUDY_DESIGN}-3

VCF_FILENAME_SHORT=${P1_JOBID}.${POPULATION}-${POP_SIZE}.5mb.vcf.gz

EXIT_STATUS=0

if [[ $EXIT_STATUS -eq "0" ]]; then
    python s2018-11-11.pipeline4.py ${VCF_FILENAME_SHORT} ${UNIQUE_ID} ${NUM_SIMS} ${NUM_PARAMS}
    EXIT_STATUS=$(( ${EXIT_STATUS} || $? ))
    if [[ $EXIT_STATUS -ne "0" ]]; then
            echo "FAIL: pipeline4.py failed"
    fi
fi

qstat -j ${JOB_ID}

END=$(date +%s)
echo end_time: $END
echo Run time was "$(expr $END - $START)" seconds

if [ $EXIT_STATUS -eq "0" ]; then
        echo "Complete"
else
        echo "Complete with Errors"
fi