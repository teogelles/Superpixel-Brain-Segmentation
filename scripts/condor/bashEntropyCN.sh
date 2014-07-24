#! /bin/bash
# foo=200
# bar=$1
# ((foo+=bar))
matlab -nojvm -nosplash -r "testEntropy('CN',$1 + 1); exit"
