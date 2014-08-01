#! /bin/bash
<<<<<<< HEAD
matlab -nojvm -nosplash -r "condor_runSLIC(1,'AD',1,250,.05+ (.01*($1)),18); exit"
=======
<<<<<<< HEAD
foo=200
bar=$1
((foo+=bar))
matlab -nojvm -nosplash -r "condor_runSLIC(1,'AD',1,250,.05+ (.01*($1)),18); exit"
=======
#foo=200
#bar=$1
#((foo+=bar))
matlab -nodisplay -r "addpath(genpath('/acmi/summer2014/tgelles1/brainseg2014/scripts')); condor_runSLIC(1,'AD',1,250,.05+(.05*($1)),1); exit"
>>>>>>> 08d03ea63b57820d148bcf48b770444709c1a63c
>>>>>>> 7a57d238b26608723433b6d424c4617f975d284e
