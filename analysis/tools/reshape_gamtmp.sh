#!/bin/bash

if test $# -gt 2
then
    echo "bash $0 input_photon_file out_directory"
    exit
fi
if test $# -lt 1
then
    echo "bash $0 input_photon_file out_directory"
    exit
fi
if test $# -eq 2
then
    ifile=$1
    odir=$2
fi
if test $# -eq 1
then
    ifile=$1
    odir="./result"
fi

mkdir -p $odir

head -n2 $ifile |tail -n1 > $odir/wava.txt
tail -n+3 | awk "NR%2==1" > $odir/header_table.txt
tail -n+3 | awk "NR%2==0" > $odir/data_table.txt
