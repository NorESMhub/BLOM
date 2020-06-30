#!/bin/bash
set -ex

EXPID=BLOM_channel_10
RUNDIR=/cluster/work/users/$USER/BLOM/run/$EXPID
CONFDIR=$HOME/BLOM/BLOM_fork/BLOM/regions_and_indices
#CONFDIR=/nird/home/$USER/MNP2/Datain_sigma_2_era
#ATMDIR=/cluster/shared/noresm/micom/NCEP
#ATMDIR=/cluster/shared/noresm/micom/ERA40
#CLIMDIR=$ATMDIR/clim
EXEDIR=$HOME/BLOM/BLOM_fork/BLOM/build
SUBDIR=$HOME/BLOM/BLOM_fork/BLOM/run
NTASKS=32
NTHREADS=1

#rm -rf $RUNDIR
mkdir -p $RUNDIR
cd $RUNDIR
cp $CONFDIR/*.nc .
cp $CONFDIR/*.dat .
#cp $ATMDIR/*.nc .
#cp $CLIMDIR/* .
cp -u $SUBDIR/limits .
cp -u $EXEDIR/blom .

  cat <<EOF >batchscript.sh
#!/bin/bash 
#
#  This script will launch micom
#
#SBATCH --account=nn1002k
#SBATCH --job-name=BLOM.$EXPID
#SBATCH --time=48:00:00
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=$NTASKS
#SBATCH --cpus-per-task=$NTHREADS

export OMP_NUM_THREADS=$NTHREADS

module purge --force
module load StdEnv
module load intel/2018a
module load netCDF-Fortran/4.4.4-intel-2018a-HDF5-1.8.19
module load PnetCDF/1.8.1-intel-2018a
module list

cd $RUNDIR
mpirun ./blom

EOF

chmod 755 batchscript.sh
sbatch batchscript.sh
