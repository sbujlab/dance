#! /bin/sh

while IFS= read -r line; do
    runnum=$line;
  #   ./dance -r $runnum -c conf/devika_bpm-regression_early_runs.conf
    ./dance -r $runnum -c conf/devika_bpm-regression_early_runs_n1Y.conf
done < $1
