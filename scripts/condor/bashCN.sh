#! /bin/bash
foo=200
bar=$1
((foo+=bar))
matlab -nojvm -nosplash -r "condor_runSLICExact('CN', $1 + 1, .1, 18); exit"
