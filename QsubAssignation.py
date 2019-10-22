# coding: utf-8
"""Python3.6"""
# compatibility: python2.7, python2.6

import time
from optparse import OptionParser
import os

sCurrentVersionScript="v1"
iTime1=time.time()
########################################################################
'''
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
# TASK_MAX=100

CONF_COMMENT="#"
CONF_STEP="="
KEYCONF_SCALL="SCALL"
KEYCONF_SPARAM="SPARAM"
KEYCONF_STASKARRAY="STASKARRAY"
KEYCONF_SMAXTASK="SMAXTASK"
KEYCONF_SMAXSIMJOB="SMAXSIMJOB"
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

########################################################################
#Function 	
def LoadConfFile(sPath):
	dDict={}
	for sNewLine in open(sPath):
		# sLine=sNewLine.strip()
		sLine=sNewLine.replace("\n","")
		sLine=sLine.replace("\\","")
		if len(sLine)==0:
			continue
		if sLine[0]==CONF_COMMENT:
			continue
		tLine=sLine.split(CONF_STEP)
		dDict[tLine[0]]=CONF_STEP.join(tLine[1:])
	return dDict	
	
def GetListFile(sDir):
	tAllFile=os.listdir(sDir)
	tFile=[]
	for sFileName in sorted(tAllFile):
		sFile=sFileName[1:] #replace 1R_/2R_ by R_
		if len(tFile)>0:
			if sFile in tFile:
				continue
		tFile.append(sFile)
	return tFile
		
def WriteBash(tList,sScriptDir,sKmerPath,sOutputPath,sDir,dCall,sConf):
	FILE=open(sOutputPath,"w")
	iSize=len(tList)
	FILE.write("#! /bin/bash\n\n")
	FILE.write(dCall[KEYCONF_SCALL]+" "+dCall[KEYCONF_SPARAM]+" "+dCall[KEYCONF_STASKARRAY]+"1-"+str(iSize)+dCall[KEYCONF_SMAXTASK]+dCall[KEYCONF_SMAXSIMJOB]+" -e "+BASHSCRIPT+".e"+dCall[KEYCONF_SPSEUDOTASKID]+" -o "+BASHSCRIPT+".o"+dCall[KEYCONF_SPSEUDOTASKID]+" "+sScriptDir+"/"+BASHSCRIPT+" "+sKmerPath+" "+sDir+" "+sScriptDir+" Demultiplexing_Ok "+sConf+"\n")
	FILE.write("""
if [ ! -d "Demultiplexing_Ok" ] ; then mkdir "Demultiplexing_Ok" ; fi
while true ; do
	if [ $(ls Demultiplexing_Ok/ | wc -l) -eq 0 ]
		then
		nbr_ok=0
	else
		nbr_ok=$(ls Demultiplexing_Ok/*_MakeAssignation.ok | wc -l)
	fi
	if [ "${nbr_ok}" -eq """+str(iSize)+""" ]
		then
		rm -r Demultiplexing_Ok
		break
	fi
	sleep 60
done\n""")
	FILE.write("echo \""+BASHSCRIPT+"\" DONE\n")
	FILE.close()
	
########################################################################
#MAIN
if __name__ == "__main__":
	dConf=LoadConfFile(sConf)
	tListFile=GetListFile(sWorkDir)
	WriteBash(tListFile,sScript,sKmerList,sOutput,sWorkDir,dConf,sConf)
	
########################################################################    
iTime2=time.time()
iDeltaTime=iTime2-iTime1
print("Script done: "+str(iDeltaTime))
