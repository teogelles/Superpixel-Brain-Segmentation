#! /bin/bash
foo=200
bar=$1
((foo+=bar))
matlab -nojvm -nosplash -r "CRFGM_test($1,???,1,true); exit"
