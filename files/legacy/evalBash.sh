#! /bin/bash
foo=8
bar=$1
((foo+=bar))
matlab -nojvm -nosplash -r "addpath(genpath('/acmi/chris13/scripts'));evalCRF($1,1); exit"
