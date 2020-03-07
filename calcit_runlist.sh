#! /bin/sh

while IFS= read -r line; do
    runnum=$line;
    ./calcit -r $runnum -c conf/sens.5376-.conf
done < $1
