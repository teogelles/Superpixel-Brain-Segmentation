#! /bin/bash
foo=200
bar=$1
((foo+=bar))
matlab -nojvm -nosplash -r "rmpath('/acmi/summer2014/agilchr1/brainseg2014/scripts/');addpath(genpath('/acmi/chris13/scripts'));sm_CRFGM_test($1,280,1); exit"
