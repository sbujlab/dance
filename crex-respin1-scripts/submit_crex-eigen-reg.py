#!/apps/python/PRO/bin/python
from subprocess import call
import subprocess
import sys,os,time

def main():
    
    _email="cameronc@jlab.org"
    _mssdir="/mss/halla/parity/raw"
    _source="/u/group/halla/parity/software/japan_offline/lagrange/crex-respin1-scripts"
    _directory="/lustre/expphy/cache/halla/parity/raw"
    _rootout="./rootfiles/"
    #_nrStart=8558
    _nrStart=5000
    _nrStop=9000
    submit=1
    useSWIF=1 #0: uses jsub 1: uses SWIF+jsub

    firstrun=9999
    lastrun=0
    _runlist=[]
    #runfile=open(_source+"/prex-runlist/simple_list/all.list","r")
    runfile=open("/u/group/halla/parity/software/japan_offline/prompt/prex-prompt/crex-runlist/simple_list/all_crex.list","r")
    for line in runfile:
        _runlist.append(int(line))
        if (firstrun >= int(line) and _nrStart <= int(line)):
            firstrun=int(line)
        if (lastrun <= int(line) and _nrStop >= int(line)):
            lastrun=int(line)
    runfile.close()
                
    _workflowID="crex-eigen-reg-parts_"+str(firstrun)+"_"+str(lastrun)

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
    f.write("  <Track name=\"one_pass\"/>\n")
    f.write("  <Name name=\""+workflowID+"\"/>\n")
    f.write("  <OS name=\"centos77\"/>\n")
    f.write("  <Memory space=\"2000\" unit=\"MB\"/>\n")
    for nr in runlist:
        if (nr < nStart or nr > nStop):
            continue
        f.write("  <Job>\n")
        f.write("    <Command><![CDATA[\n")
        f.write("    cd "+source+"\n")
        f.write("    cd ../\n")
        f.write("    source /site/12gev_phys/softenv.csh 2.4\n") # Tao updated the most recent compilation to 2.4
        f.write("    "+source+"/crex-eigen-reg-parts.sh "+str(nr)+"; \n")
        f.write("    ]]></Command>\n")
        f.write("    <Stdout dest=\""+source+"/LogFiles/crex-eigen-reg_respin1log"+"_%04d"%(nr)+".out\"/>\n")
        f.write("    <Stderr dest=\""+source+"/LogFiles/crex-eigen-reg_respin1log"+"_%04d"%(nr)+".err\"/>\n")
        f.write("  </Job>\n\n")

    f.write("</Request>\n")
    f.close()
    return 0
                    
if __name__ == '__main__':
    main()
