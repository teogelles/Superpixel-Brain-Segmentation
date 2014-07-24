#! /bin/bash
foo=200
bar=$1
((foo+=bar))
matlab -nojvm -nosplash -r "condor_runSLIC($1 + 1, 'CN'); exit"
