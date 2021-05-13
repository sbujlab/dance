#!/bin/bash
filelist=""
shopt -s extglob
while IFS= read -r line; do
    all_files=$(ls -1 ./dit-coeffs/prexPrompt_ditcoeffs_$line.000.root);
    for file in $all_files
    do
	echo $file;
	filelist+=" "$file;
    done

done < ./prex-runlist/simple_list/slug$1.list
echo $filelist;

hadd -f ./ditcoeffs_by_slug/slug$1.root $filelist;

shopt -u extglob
