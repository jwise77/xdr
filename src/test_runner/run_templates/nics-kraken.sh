#!/bin/bash
#PBS -N ${TEST_NAME}
#PBS -j oe
#PBS -l walltime=${RUN_TIME},size=${N_PROCS}
#PBS -A ${ACCOUNT}

cd $PBS_O_WORKDIR

aprun -n ${N_PROCS} ${ENZO_EXE} -d ${PAR_FILE} >& estd.out;