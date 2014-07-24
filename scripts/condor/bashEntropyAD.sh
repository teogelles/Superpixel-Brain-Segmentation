#! /bin/bash
# foo=200
# bar=$1
# ((foo+=bar))
matlab -nojvm -nosplash -r "testEntropy('AD',$1 + 1); exit"
