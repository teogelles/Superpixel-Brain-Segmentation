#! /bin/bash
foo=200
bar=$2
((foo+=bar))
matlab -nojvm -nosplash -r "condor_runSLICExact($1,$2 + 1); exit"
