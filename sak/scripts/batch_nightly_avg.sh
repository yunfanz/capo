#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -l h_vmem=2G
#$ -j y
#$ -N nightly_xtalk_xx
#$ -o grid_output

CAL=psa6240_v003
source /usr/global/paper/CanopyVirtualEnvs/PAPER_Polar/bin/activate
POL='xx'

JDS=`python list_jds.py $*`
MYJDS=`pull_args.py ${JDS}`
echo working on JDS: $MYJDS
for JD in $MYJDS
do
    MYFILES=`python -c "print ' '.join([l for l in '''$*'''.split() if l.find('${JD}')>0])"`
    echo Working on $MYFILES
    python nightly_avg.py $MYFILES
done
