#!/bin/sh

root -b -q 'rootScripts/ProbeCycle_prex.C';
root -b -q 'rootScripts/AverageSensitivity.C';
root -b -q 'rootScripts/AverageSlope.C';
root -b -q 'rootScripts/SolveMergedCycles.C';
root -b -q 'rootScripts/ResidualSensByCyle.C';
root -b -q 'rootScripts/ResidualSensByRun.C';

