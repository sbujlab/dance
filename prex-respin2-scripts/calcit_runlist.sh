#! /bin/sh

while IFS= read -r line; do
    runnum=$line;
    ./auto_calcit.sh $runnum;
done < $1
