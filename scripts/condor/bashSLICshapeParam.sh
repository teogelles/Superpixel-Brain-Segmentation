#! /bin/bash
#foo=200
#bar=$1
#((foo+=bar))
matlab -nodisplay -r "addpath(genpath('/acmi/summer2014/tgelles1/brainseg2014/scripts')); condor_runSLIC(1,'AD',1,250,.05+(.05*($1)),1); exit"
