#! /bin/bash
foo=8
bar=$1
((foo+=bar))
matlab -nojvm -nosplash -r "addpath(genpath('/acmi/chris13/scripts'));CRFGM_ADNI('/scratch/cmagnan1/testing7res1num1/paramsIter528.mat',$1,1,1); exit"
