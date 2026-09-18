#!/bin/tcsh -f
setenv STAGE2_MANIFEST $1
setenv STAGE2_TAG $2
root4star -b -q Stage2_BatchShard.C+

