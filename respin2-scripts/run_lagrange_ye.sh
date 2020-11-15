#!/bin/bash
for i in `seq $1 $2`; 
do
  echo $i
  ./run_lagrange_slug_ye.sh  ./prex-runlist/simple_list/slug$i.list 
done
