#! /bin/bash
foo=8
bar=$1
((foo+=bar))
matlab -nojvm -nosplash -r "addpath(genpath('/acmi/chris13/scripts'));testv1CRF($1,1); exit"
