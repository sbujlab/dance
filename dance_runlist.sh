#! /bin/sh

while IFS= read -r line; do
    runnum=$line;
    # ./dance -r $runnum -c conf/loadtest.conf
    ./lagrange.sh $runnum;
done < $1
