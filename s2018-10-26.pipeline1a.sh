#$ -S /bin/sh
#$ -cwd
#$ -r y
#$ -o /netapp/home/dmctong/imp_rvat/50k_scripts/2019-01/
#$ -e /netapp/home/dmctong/imp_rvat/50k_scripts/2019-01/
#$ -l mem_free=4G
#$ -l h_rt=336:00:00
#$ -t 1
#$ -l arch=linux-x64

echo SGE_TASK_ID: $SGE_TASK_ID
echo JOB_ID: $JOB_ID

i=$(expr $SGE_TASK_ID - 1)
echo i: $i

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/gsl-2.4/lib

echo $LD_LIBRARY_PATH

START=$(date +%s)
echo start_time: $START

P1_JOBID=$1
POPULATION=$2
POP_SIZE=$3

python /netapp/home/dmctong/imp_rvat/50k_scripts/s2018-10-26.pipeline1a.py ${P1_JOBID} ${POPULATION} ${POP_SIZE}

END=$(date +%s)
echo end_time: $END
echo Run time was "$(expr $END - $START)" seconds