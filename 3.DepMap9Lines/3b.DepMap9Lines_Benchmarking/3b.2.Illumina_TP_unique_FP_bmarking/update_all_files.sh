#!/bin/bash

set -e

cwd=`pwd`

for dir in "__illum_TP_uniq_FP.starF" "__illum_TP_uniq_FP.arriba,starF" "__illum_TP_uniq_FP.either" "__illum_TP_uniq_FP.arriba"; do
    echo updating $dir
    cd $dir/data
    ./update_files.sh
    cd $cwd
done
