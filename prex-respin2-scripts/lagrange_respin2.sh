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

    conf_tr="lagrange-respin2-trunc.conf" 
    conf1="lagrange-respin2-eigen.conf"
    conf2="lagrange-respin2-forced.conf"
    conf3="lagrange-respin2.conf"
    ./dance \
    	-f $rootfile \
    	-c ./conf/$conf_tr ; 
    # ./dance \
    # 	-f $rootfile \
    # 	-c ./conf/$conf1 ; 
    # ./dance \
    # 	-f $rootfile \
    # 	-c ./conf/$conf2 ; 
    # ./dance \
    # 	-f $rootfile \
    # 	-c ./conf/$conf3 ; 
done

