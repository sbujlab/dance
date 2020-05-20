#!/bin/bash
filelist=""
while IFS= read -r line; do
    filename=" ./dit-coeffs/prexPrompt_ditcoeffs_"$line".*.root";
    echo $filename;
    if [ -f $filename ]; then
	filelist+=$filename;
    fi
done < ./prex-runlist/simple_list/slug$1.list
echo $filelist;

hadd -f ./ditcoeffs_by_slug/slug$1.root $filelist;
