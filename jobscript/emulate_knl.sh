#! /bin/bash
#SBATCH -J KNL-emu
#SBATCH -p mixed
#SBATCH -N 1
#SBATCH -t 12:00:00
#SBATCH --ntasks-per-node=1
#SBATCH -o joblog/emulate_knl.explicit.log

cd $SLURM_SUBMIT_DIR

module purge
module load intel intelmpi mkl

export OMP_STACKSIZE="2M"
export OMP_SCHEDULE="static"
/work/TCAPROFC/hirokawa/sde/sde-external-7.41.0-2016-03-03-lin/sde -knl -- ./bin/ARTED_sc.mic < data/input_sc_Si.single.dat

