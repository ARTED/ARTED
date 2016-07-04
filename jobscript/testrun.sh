#! /bin/bash
#SBATCH -J armic-native
#SBATCH -p debug
#SBATCH -N 1
#SBATCH --ntasks-per-node=2
#SBATCH -t 00:30:00
#SBATCH -o joblog/testrun.log

cd $SLURM_SUBMIT_DIR

module purge
module load intel intelmpi mkl
#module load openmpi lapack

export OMP_STACKSIZE="2M"
export OMP_SCHEDULE="static"
export MIC_KMP_AFFINITY="scatter,granularity=fine"

TARGET_MODE=sc
INPUT_FILE=./data/input_${TARGET_MODE}_Si.dat
#INPUT_FILE=./data/input_${TARGET_MODE}_Si.single.dat
#INPUT_FILE=./data/input_${TARGET_MODE}.dat

#export MIC_NPN=1

#export MIC_OMP_PLACES="4:240"
#export MIC_OMP_PROC_BIND=spread,close
#export MIC_KMP_AFFINITY="explicit,proclist=[4-243]"
#export MIC_OMP_NESTED=1
#export MIC_OMP_THREAD_LIMIT=240
#export MIC_OMP_MAX_ACTIVE_LEVELS=2
#export MIC_OMP_NUM_THREADS=240
#export MIC_MKL_DYNAMIC=0
#export MIC_ARTED_OUTER_THREADS=8
#export MIC_ARTED_INNER_THREADS=30

#mpirun --map-by numa \
#       ./script/env_wrapper ./bin/ARTED_${TARGET_MODE} \
#       < ${INPUT_FILE}

./script/mpirun.symm -c ./bin/ARTED_${TARGET_MODE} < ${INPUT_FILE}
#./script/mpirun.symm -c ./bin/ARTED_${TARGET_MODE} < ${INPUT_FILE}
#./script/mpirun.symm -c ./bin/ARTED_${TARGET_MODE} < ${INPUT_FILE}

#export MIC_OMP_NUM_THREADS=120
#./script/mpirun.symm -m ./bin/ARTED_${TARGET_MODE} < ${INPUT_FILE}
#./script/mpirun.symm -m ./bin/ARTED_${TARGET_MODE} < ${INPUT_FILE}
#./script/mpirun.symm -m ./bin/ARTED_${TARGET_MODE} < ${INPUT_FILE}
#
#export MIC_OMP_NUM_THREADS=180
#./script/mpirun.symm -m ./bin/ARTED_${TARGET_MODE} < ${INPUT_FILE}
#./script/mpirun.symm -m ./bin/ARTED_${TARGET_MODE} < ${INPUT_FILE}
#./script/mpirun.symm -m ./bin/ARTED_${TARGET_MODE} < ${INPUT_FILE}
#
export MIC_OMP_NUM_THREADS=240
./script/mpirun.symm -m ./bin/ARTED_${TARGET_MODE} < ${INPUT_FILE}
#./script/mpirun.symm -m ./bin/ARTED_${TARGET_MODE} < ${INPUT_FILE}
#./script/mpirun.symm -m ./bin/ARTED_${TARGET_MODE} < ${INPUT_FILE}
