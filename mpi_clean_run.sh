#!/bin/bash

NUM_PROCS=4

make clean
make -j 4
mpirun -np $NUM_PROCS rtdVal > run_log