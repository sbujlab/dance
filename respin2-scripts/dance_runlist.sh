#! /bin/sh

while IFS= read -r line; do
    runnum=$line;
    ./updownDD_respin2.sh $runnum;
done < $1
