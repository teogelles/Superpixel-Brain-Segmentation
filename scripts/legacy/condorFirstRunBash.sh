#! /bin/bash
foo=200
bar=$1
((foo+=bar))
matlab -nojvm -nosplash -r "CRF_test($1,280,1,true); exit"
