#!/bin/ksh
set -ex

EXPID=C05
RUNDIR=/work/users/matsbn/MNP4/run/$EXPID
CONFDIR=/home/nersc/matsbn/MNP4/Datain_sigma_2
ATMDIR=/work/shared/nersc/gcr/NCEP
CLIMDIR=$ATMDIR/clim
EXEDIR=/home/nersc/matsbn/MNP4/Micom/$EXPID/build
PARTDIR=/home/nersc/matsbn/MNP4/Partit
SUBDIR=/home/nersc/matsbn/MNP4/Micom/$EXPID/run
NPROC=231

#rm -rf $RUNDIR
mkdir -p $RUNDIR
lfs setstripe $RUNDIR 0 -1 -1
cd $RUNDIR
cp -u $CONFDIR/* .
cp -u $ATMDIR/*.nc .
cp -u $CLIMDIR/* .
cp -u $SUBDIR/limits .
cp $PARTDIR/patch.input.$NPROC patch.input
cp -u $EXEDIR/micom .

  cat <<EOF >qsub.sh
#!/bin/ksh 
#
#  This script will launch micom
#
#PBS -A nn2345k
#PBS -N micom.$EXPID
#PBS -l mppwidth=$NPROC
#PBS -l mppnppn=4
#PBS -l walltime=06:00:00
#PBS -j oe

cd $RUNDIR
aprun -n $NPROC -N 4 ./micom

EOF

chmod 755 qsub.sh
qsub qsub.sh
