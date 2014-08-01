#! /bin/bash
foo=200
bar=$1
((foo+=bar))
../condor_tests/run_condor_runSLICExact.sh /usr/local/MATLAB/MATLAB_Compiler_Runtime/v83 MCI $1+1 0.1 18
