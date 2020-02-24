# coding: utf-8
"""Python3.6"""
# compatibility: python2.7, python2.6

import time
from optparse import OptionParser
import os

sCurrentVersionScript="v3"
iTime1=time.time()
########################################################################
'''
V4-2020/02/11
Adapt to multiplate analysis

V3-2020/01/21
Adapt to array wioth limited value
V2-2019/10/30
Work with index on base file instead multiple subfile (decrease memory usage)
V1-2019/07/01
Prepare Qsub command list for Assignation task

python QsubAssignation.py -s SCRIPT_ASSIGNATION -k KMER_LIST -d DEMULTIPLEXING_DIR -o OUTPUT
SCRIPT_ASSIGNATION: Script that make the assignation
KMER_LIST: Tabulated file that contains Kmer to sampleId
DEMULTIPLEXING_DIR: Work directory for the assignation script
OUTPUT: Script ready to use
'''
########################################################################
#CONSTANT
BASHSCRIPT="LaunchAssignation.sh"

DEFAULT_SEQ_BY_TASK=1000000 #Must be similar in MakeAssignation.py

CONF_COMMENT="#"
CONF_STEP="="
KEYCONF_SCALL="SCALL"
KEYCONF_SPARAM="SPARAM"
KEYCONF_STASKARRAY="STASKARRAY"
KEYCONF_SMAXTASK="SMAXTASK"
KEYCONF_SMAXSIMJOB="SMAXSIMJOB"
KEYCONF_SMAXARRAYSIZE="SMAXARRAYSIZE"
KEYCONF_STASKID="STASKID"
KEYCONF_SPSEUDOTASKID="SPSEUDOTASKID"
########################################################################
#Options
parser = OptionParser()
parser.add_option("-s","--script", dest="script")
parser.add_option("-k","--kmerlist", dest="kmerlist")
parser.add_option("-d","--workdir", dest="workdir")
parser.add_option("-o","--output", dest="output")
parser.add_option("-c","--conffile", dest="conf")
parser.add_option("-q","--quantity", dest="quantity")
parser.add_option("-a","--argfile", dest="argfile")
parser.add_option("-p","--pid", dest="pid")
parser.add_option("-v","--varname", dest="varname")

(options, args) = parser.parse_args()

sScript=options.script
if not sScript:
	exit("Error : no script -s defined, process broken")
	if sScript[-1]=="/":
		sScript=sScript[:-1]

sKmerList=options.kmerlist
if not sKmerList:
	exit("Error : no kmerlist -k defined, process broken")

sWorkDir=options.workdir
if not sWorkDir:
	exit("Error : no workdir -d defined, process broken")

sOutput=options.output
if not sOutput:
	exit("Error : no output -o defined, process broken")

sConf=options.conf
if not sConf:
	exit("Error : no conf -c defined, process broken")

sArg=options.argfile
if not sArg:
	exit("Error : no argfile -q defined, process broken")

sQuantity=options.quantity
if not sQuantity:
	exit("Error : no quantity -q defined, process broken")
try:
	iQuantity=int(sQuantity)
	# if iQuantity%SEQ_BY_TASK==0:
		# iQuantity=iQuantity/SEQ_BY_TASK
	# else:
		# iQuantity=int(round(iQuantity/SEQ_BY_TASK,0))+1
except ValueError:
	exit("Error : quantity -q must be an integer, process broken")

sPID=options.pid
if not sPID:
	exit("Error : no pid -p defined, process broken")

sVarName=options.varname
if not sVarName:
	exit("Error : no varname -v defined, process broken")

########################################################################
#Function 	
def LoadConfFile(sPath):
	dDict={}
	for sNewLine in open(sPath):
		# sLine=sNewLine.strip()
		sLine=sNewLine.replace("\n","")
		sLine=sLine.replace("\\","")
		sLine=sLine.replace('"','')
		if len(sLine)==0:
			continue
		if sLine[0]==CONF_COMMENT:
			continue
		tLine=sLine.split(CONF_STEP)
		dDict[tLine[0]]=CONF_STEP.join(tLine[1:])
	return dDict	
		
def WriteBash(sArg,iSize,sScriptDir,sKmerPath,sOutputPath,sDir,dCall,sConf,sPID,iNumberSeq,sPlateId):
	sLogDir=sPID+"_"+sPlateId+"_log_LaunchAssignation"
	FILE=open(sOutputPath,"w")
	FILE.write("#! /bin/bash\n\n")
	FILE.write("mkdir "+sLogDir+"\n")
	FILE.write("if [ ! -d Demultiplexing"+sPlateId+"_Ok ] ; then mkdir Demultiplexing"+sPlateId+"_Ok ; fi")
	FILE.write(dCall[KEYCONF_SCALL]+" "+dCall[KEYCONF_SPARAM]+" "+dCall[KEYCONF_STASKARRAY]+"1-"+str(iSize)+dCall[KEYCONF_SMAXTASK]+dCall[KEYCONF_SMAXSIMJOB]+" -e "+sLogDir+"/"+BASHSCRIPT.replace(".sh","")+".e"+dCall[KEYCONF_SPSEUDOTASKID]+" -o "+sLogDir+"/"+BASHSCRIPT.replace(".sh","")+".o"+dCall[KEYCONF_SPSEUDOTASKID]+" "+sScriptDir+"/"+BASHSCRIPT+" "+str(iNumberSeq)+" "+sKmerPath+" "+sDir+" "+sScriptDir+" Demultiplexing"+sPlateId+"_Ok "+sConf+" "+sArg+" "+sPlateId+"\n")
	# FILE.write("""
# if [ ! -d "Demultiplexing"""+sPlateId+"""_Ok" ] ; then mkdir "Demultiplexing"""+sPlateId+"""_Ok" ; fi
# while true ; do
	# if [ $(ls Demultiplexing"""+sPlateId+"""_Ok/ | wc -l) -eq 0 ]
		# then
		# nbr_ok=0
	# else
		# nbr_ok=$(ls Demultiplexing"""+sPlateId+"""_Ok/*_MakeAssignation.ok | wc -l)
	# fi
	# if [ "${nbr_ok}" -eq """+str(iSize)+""" ]
		# then
		# rm -r Demultiplexing"""+sPlateId+"""_Ok
		# break
	# fi
	# sleep 60
# done\n""")
	FILE.write("rm -r Demultiplexing"+sPlateId+"_Ok\n")
	FILE.write("echo \""+BASHSCRIPT+"\" DONE\n")
	FILE.close()

def GetJobByTask(iSeq,iTask):
	if iTask==0:
		return DEFAULT_SEQ_BY_TASK,int(round(iSeq/DEFAULT_SEQ_BY_TASK,0))+1
	if iSeq%iTask==0:
		iResult=iSeq/iTask
	else:
		iResult=int(round(iSeq/iTask,0))+1
	return iResult,iTask

########################################################################
#MAIN
if __name__ == "__main__":
	dConf=LoadConfFile(sConf)
	iJobByTask,iTask=GetJobByTask(iQuantity,int(dConf[KEYCONF_SMAXARRAYSIZE]))
	WriteBash(sArg,iTask,sScript,sKmerList,sOutput,sWorkDir,dConf,sConf,sPID,iJobByTask,sVarName)
	
########################################################################    
iTime2=time.time()
iDeltaTime=iTime2-iTime1
print("Script done: "+str(iDeltaTime))
