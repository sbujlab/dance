#! /bin/sh
runnum=$1;

if [ -z "$runnum" ] 
then
    echo "Run Number is empty";
    exit 1;
fi    

level="Prompt"			
shopt -s extglob
rootfile_list=$(ls -1 $QW_ROOTFILES/prex$level\_pass1_$runnum.!(*jlab*).root);
shopt -u extglob

for rootfile  in $rootfile_list
do
    trim=${rootfile%.root}
    run_dot_seg=${trim#*pass1_}
    run_num=${run_dot_seg%.*}
    run_seg=${run_dot_seg/./_}

    sensConf="sens.conf"
    if [ $(($run_num)) -ge 3130 ]
    then
      sensConf="sens.3130-.conf"
    fi
    if [ $(($run_num)) -ge 3390 ]
    then
      sensConf="sens.3390-.conf"
    fi
    if [ $(($run_num)) -ge 3803 ]
    then
      sensConf="sens.3803-.conf"
    fi
    if [ $(($run_num)) -ge 5376 ]
    then
      sensConf="sens.5376-.conf"
    fi
    ./calcit \
    	-f $rootfile \
    	-c ./conf/$sensConf ; 

done

