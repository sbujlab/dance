#!/apps/python/PRO/bin/python
from subprocess import call
import subprocess
import sys,os,time

def main():
    
    _email="tao@jlab.org"
    _mssdir="/mss/halla/parity/raw"
    _source="/u/home/tao/work/respin2-sandbox/dance_respin2/"
    _directory="/lustre/expphy/cache/halla/parity/raw"
    _rootout="/u/home/tao/work/respin2-sandbox/dance_respin2/rootfiles/"
    _nrStart=2000
    _nrStop=6000
    submit=1
    useSWIF=1 #0: uses jsub 1: uses SWIF+jsub

    firstrun=9999
    lastrun=0
    _runlist=[]
    # runfile=open(_source+"/prex-runlist/simple_list/extra_march9.list","r")
    # runfile=open(_source+"/prex-runlist/good_list/test.list","r")
    #runfile=open(_source+"/prex-runlist/simple_list/extra_feb28.list","r")
    runfile=open(_source+"/prex-runlist/simple_list/all.list","r")
    # runfile=open(_source+"/prex-runlist/simple_list/slug1.list","r")
    for line in runfile:
        _runlist.append(int(line))
        if (firstrun >= int(line) and _nrStart <= int(line)):
            firstrun=int(line)
        if (lastrun <= int(line) and _nrStop >= int(line)):
            lastrun=int(line)
    runfile.close()
                
    _workflowID="Calcit_"+str(firstrun)+"_"+str(lastrun)

    createXMLfile(_mssdir,_source,_rootout,_nrStart,_nrStop,_email,_workflowID,_runlist)

    if submit==1:
        if useSWIF==1:
            print "Submitting lagrange prompt analysis for runs "+str(firstrun)+" to "+str(lastrun)+" using designated SWIF workflow "+str(_workflowID)
            call(["swif","add-jsub","-workflow",str(_workflowID),"-create","-script",_source+"/"+_workflowID+".xml"])
        elif useSWIF==0:
            print "submitting position sampled with id between ",firstrun,lastrun
            call(["jsub","-xml",_source+"/"+_workflowID+".xml"])
        else:
            print "NOT submitting position sampled with id between ",firstrun,lastrun
            
    print "I am all done"


def createXMLfile(mssdir,source,rootout,nStart,nStop,email,workflowID,runlist):

    f=open(source+"/"+workflowID+".xml","w")
    f.write("<Request>\n")
    f.write("  <Email email=\""+email+"\" request=\"false\" job=\"true\"/>\n")
    f.write("  <Project name=\"prex\"/>\n")
    #    f.write("  <Track name=\"debug\"/>\n")
    f.write("  <Track name=\"one_pass\"/>\n")
    f.write("  <Name name=\""+workflowID+"\"/>\n")
    f.write("  <OS name=\"centos7\"/>\n")
    f.write("  <Memory space=\"2000\" unit=\"MB\"/>\n")

    # for nr in range(nStart,nStop+1): # repeat for nr jobs
    # print "    <Stdout dest=\""+source+"/LogFiles/postpan_ifarmlog"+"_%04d"%(runlist[0])+".out\"/>\n"
    # print "    <Stderr dest=\""+source+"/LogFiles/postpan_ifarmlog"+"_%04d"%(runlist[0])+".err\"/>\n"
    for nr in runlist:
        if (nr < nStart or nr > nStop):
            continue
        f.write("  <Job>\n")
        cmd = "ls /mss/halla/parity/raw/parity_ALL_"+str(nr)+".dat.* | sort -t'.' -n -k3 | tail -n 1 | awk -F. '{print $3}'"
        ps = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        output = ps.communicate()[0]
        # for partfile in range(0,int(output)+1):
        #     datName="mss:"+mssdir+"/"+"parity_ALL"+"_%04d"%(nr)+".dat."+str(partfile)
        # f.write("    <Input src=\""+datName+"\" dest=\"parity_ALL"+"_%04d"%(nr)+".dat."+str(partfile)+"\"/>\n")
        f.write("    <Command><![CDATA[\n")
        f.write("    setenv QW_ROOTFILES /lustre/expphy/volatile/halla/parity/japanOutput/ \n");
        f.write("    cd "+source+"\n")
        f.write("    "+source+"/auto_calcit.sh "+str(nr)+"; \n")
        f.write("    ]]></Command>\n")
        f.write("    <Stdout dest=\""+source+"/LogFiles/calcit_ifarmlog"+"_%04d"%(nr)+".out\"/>\n")
        f.write("    <Stderr dest=\""+source+"/LogFiles/calcit_ifarmlog"+"_%04d"%(nr)+".err\"/>\n")
        f.write("  </Job>\n\n")

    f.write("</Request>\n")
    f.close()
    return 0
                    
if __name__ == '__main__':
    main()
