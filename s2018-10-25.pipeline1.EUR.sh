#$ -S /bin/sh
#$ -cwd
#$ -r y
#$ -o /netapp/home/dmctong/imp_rvat/50k_scripts/2019-01/
#$ -e /netapp/home/dmctong/imp_rvat/50k_scripts/2019-01/
#$ -l mem_free=3G
#$ -l h_rt=336:00:00
#$ -t 1

echo SGE_TASK_ID: $SGE_TASK_ID
echo JOB_ID: $JOB_ID

i=$(expr $SGE_TASK_ID - 1)
echo i: $i

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/gsl-2.4/lib
echo $LD_LIBRARY_PATH

START=$(date +%s)
echo start_time: $START

POP_SIZE=50000
echo ${POP_SIZE}

EXIT_STATUS=0

if [[ $EXIT_STATUS -eq "0" ]]; then
    TEMP_DIR=`mktemp -d -p /scratch/`
    EXIT_STATUS=$(( ${EXIT_STATUS} || $? ))
    if [[ $EXIT_STATUS -ne "0" ]]; then
            echo "FAIL: Making working directory failed"
    fi
fi

# python /netapp/home/dmctong/imp_rvat/scripts/s2018-09-13.pipeline1.EUR.py $POP_SIZE $JOB_ID /wynton/scratch/dmctong/imp_rvat/data/input-gzip/${JOB_ID}.EUR-${POP_SIZE}.vcf.gz
python /netapp/home/dmctong/imp_rvat/50k_scripts/s2018-10-25.pipeline1.EUR.py $POP_SIZE $JOB_ID ${JOB_ID}.EUR-${POP_SIZE}.vcf.gz $TEMP_DIR

# rsync back to here:

mkdir /wynton/scratch/dmctong/imp_rvat/data/input-gzip/ # just in case
mkdir /wynton/scratch/dmctong/imp_rvat/data/input-gzip/EUR/ # just in case
if [[ $EXIT_STATUS -eq "0" ]]; then
    rsync -rauz ${TEMP_DIR}/* /wynton/scratch/dmctong/imp_rvat/data/input-gzip/EUR/
    EXIT_STATUS=$(( ${EXIT_STATUS} || $? ))
    if [[ $EXIT_STATUS -ne "0" ]]; then
            echo "FAIL: Return transfer failed"
    fi
fi

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