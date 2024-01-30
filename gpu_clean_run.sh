#!/bin/bash

make clean
make -j 6
./rtdVal | tee run_log
