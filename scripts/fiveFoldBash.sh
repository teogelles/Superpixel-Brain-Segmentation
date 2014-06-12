#! /bin/bash
matlab -nojvm -nosplash -r "addpath(genpath('/home/cmagnan1/MRISegmentation/UGM'));addpath(genpath('/acmi/fmri/spm8'));addpath(genpath('/acmi/fmri/UGM/CRFcell'));addpath(genpath('/acmi/chris13/scripts'));CRFGM2($1); exit"
