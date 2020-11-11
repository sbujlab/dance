#!/bin/bash
filelist=""
shopt -s extglob
while IFS= read -r line; do
    ./bcm_respin2.sh $line;
    all_files=$(ls -1 ./rootfiles/charge_sensitivity/prexRespin2_bcm_eigen_$line.000.root);
    
    for file in $all_files
    do
	echo $file;
	filelist+=" "$file;
    done
done < ./prex-runlist/simple_list/slug$1.list

echo $filelist;

hadd -f ./treeMergeOutput/Merged_BCM_slug$1.root $filelist;

shopt -u extglob
