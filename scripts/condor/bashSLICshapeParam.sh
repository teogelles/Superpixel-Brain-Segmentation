#! /bin/bash
foo=200
bar=$1
((foo+=bar))
matlab -nojvm -nosplash -r "condor_runSLIC(1,'AD',1,250,.05+ (.01*($1)),18); exit"
