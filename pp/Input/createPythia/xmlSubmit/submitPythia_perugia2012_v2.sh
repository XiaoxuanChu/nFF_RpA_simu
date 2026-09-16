#!/bin/bash

ipt=0
#OUTDIR=/star/data05/scratch/tinglin/simulation/pythia_perugia0/pp200
OUTDIR=/star/data05/scratch/tinglin/simulation/pythia_perugia2012/pp200
#OUTDIR=/star/data05/scratch/tinglin/simulation/pythia_tuneA/pp200
        if [ -e $OUTDIR ];then
                rm -r $OUTDIR
        fi
        mkdir -p $OUTDIR
        mkdir -p $OUTDIR/log
        mkdir -p $OUTDIR/err
for ptbin in pt2_3 pt3_4 pt4_5 pt5_7 pt7_9 pt9_11 pt11_15 pt15_25 pt25_35 pt35_-1; do
#for ptbin in pt9_11;do
    ptmin=`echo $ptbin|sed 's/pt\(.*\)_\(.*\)/\1/'`
    ptmax=`echo $ptbin|sed 's/pt\(.*\)_\(.*\)/\2/'`
    echo ptmin=$ptmin
    echo ptmax=$ptmax
    ipt=`awk -v ptbin=$ptbin '$1 == ptbin { print $2 }' nevents.txt`
    iters=`awk -v ptbin=$ptbin '$1 == ptbin { print $3 }' nevents.txt`
    echo id=$ipt
    echo iterations=$iters

    for ((cyclenumber=1; cyclenumber <= $iters; cyclenumber++)) ; do

        outfile=${ptbin}_${cyclenumber}_pp200_perugia2012
#        outfile=${ptbin}_${cyclenumber}_pp200_perugia0
#        outfile=${ptbin}_${cyclenumber}_pp200_tuneA
        run=$((ipt*10000+cyclenumber))
        echo Submit $outfile 100000 $run

        star-submit-template -template RunPythia_perugia2012_v0.xml -entities outdir=${OUTDIR},outfile=${outfile},run=${run},nevents=100000,ptmin=${ptmin},ptmax=${ptmax}
    done

done
