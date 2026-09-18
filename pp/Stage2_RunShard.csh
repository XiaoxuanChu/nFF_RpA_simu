#!/bin/tcsh -f
# Arguments: manifest-path shard-tag.  Runs inside one scheduler scratch area.
set manifest = "$1"
set tag = "$2"
root4star -b -q "Stage2_pp_Baseline.C+(\"$manifest\",\"partial_${tag}.root\",\"partial_${tag}.txt\",0)"

