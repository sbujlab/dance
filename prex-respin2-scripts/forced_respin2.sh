#! /bin/sh
runnum=$1;

if [ -z "$runnum" ] 
then
    echo "Run Number is empty";
    exit 1;
fi    

level="Prompt"
shopt -s extglob
# find split file
rootfile_list=$(ls -1 /lustre/expphy/volatile/halla/parity/prex-respin2/japanOutput/prex$level\_pass2_$runnum.000.root);
shopt -u extglob

for rootfile  in $rootfile_list
do
    trim=${rootfile%.root}
    run_dot_seg=${trim#*pass2_}
    run_num=${run_dot_seg%.*}
    run_seg=${run_dot_seg/./_}

    echo $(($run_num))
    echo $run_num

    conf="lagrange-respin2-forced.conf"
    ./dance \
    	-f $rootfile \
    	-c ./conf/$conf ; 
done

