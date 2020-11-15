#! /bin/sh

while IFS= read -r line; do
    runnum=$line;
./dance -r $runnum -c /w/halla-scifs17exp/parity/disk1/yetian/prex/prex_lagrange/conf/lagrange-respin2.conf 

done < $1
