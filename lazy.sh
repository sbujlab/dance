#!/bin/sh

slug_id=$1;
if [ -z "$slug_id" ]
then
    echo "Slug Number is empty";
    exit 1;
fi
./merge_by_slug.sh $1;
root -b -q 'rootScripts/ProbeCycle_prex.C('$1')';
root -b -q 'rootScripts/AverageSensitivity.C('$1')';
root -b -q 'rootScripts/AverageSlope.C('$1')';
root -b -q 'rootScripts/SolveMergedCycles.C('$1')';

