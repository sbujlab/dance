#! /bin/sh

while IFS= read -r line; do
    runnum=$line;
    ./dance -r $runnum -c conf/loadtest.conf
done < $1
