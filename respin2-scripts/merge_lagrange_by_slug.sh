#!/bin/bash
filelist=""
shopt -s extglob
while IFS= read -r line; do
    all_files=$(ls -1 ./rootfiles/prexRespin2_lagrange_trunc_$line.000.root);
    
    for file in $all_files
    do
	echo $file;
	filelist+=" "$file;
    done
done < ./prex-runlist/simple_list/slug$1.list

echo $filelist;

hadd -f ./treeMergeOutput/MergedLagrange_trunc_slug$1.root $filelist;

shopt -u extglob
