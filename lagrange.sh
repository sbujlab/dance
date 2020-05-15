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
rootfile_list=$(ls -1 /lustre/expphy/volatile/halla/parity/japanOutput/prex$level\_pass1_$runnum.!(*farm*).root);
shopt -u extglob

for rootfile  in $rootfile_list
do
    trim=${rootfile%.root}
    run_dot_seg=${trim#*pass1_}
    run_num=${run_dot_seg%.*}
    run_seg=${run_dot_seg/./_}

    echo $(($run_num))
    echo $run_num

    conf="lagrange.conf"
    if [ $(($run_num)) -ge 3390 ]
    then
      conf="lagrange.3390-.conf"
    fi
    if [ $(($run_num)) -ge 4735 ]
    then
      conf="lagrange.4735-.conf"
    fi

    ./dance \
    	-f $rootfile \
    	-c ./conf/$conf ; 
done

