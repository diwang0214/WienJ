#!/bin/sh

if [ ! -n "$1" ]; then
#  echo "IS NULL"
#  mkdir szJ
  sz="sz"
else
#  echo "NOT NULL"
#  mkdir $1
  sz=$1   
fi


case=`basename $(pwd)`
mkdir $sz
mkdir $sz/data
cp $case."struct" ./$sz/dat.struct
cp $case."klist" ./$sz/dat.klist
cp $case."almblmup" ./$sz/dat.almblmup
cp $case."almblmdn" ./$sz/dat.almblmdn
cp $case."r2v" ./$sz/dat.r2v
cp $case."r2vdn" ./$sz/dat.r2vdn
cp $case."radwfup" ./$sz/dat.radwfup
cp $case."radwfdn" ./$sz/dat.radwfdn
cp $case."vorbup" ./$sz/dat.vorbup
cp $case."vorbdn" ./$sz/dat.vorbdn
cp $WIENJROOT/dat.in ./$sz/dat.in
cp $WIENJROOT/sub_wien2k.lsf ./$sz/sub_wien2k.lsf
