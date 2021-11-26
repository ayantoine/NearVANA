# coding: utf-8
"""Python3.6"""

import time
from optparse import OptionParser
import os
import re
import subprocess

sCurrentVersionScript="v1"
iTime1=time.time()
########################################################################
'''
V1-2019/09/10

python DrawGaphDensity.py -a ASSEMBLY_STAT -i IDENTIFICATION_STAT -f FOLDER_STAT -o OUTPUT -p PREFIX

ASSEMBLY_STAT: NearVANA assembly stat tsv
IDENTIFCATION_STAT: NearVANA identification stat tsv
FOLDER_STAT: NearVANA folder stat
OUTPUT: output files
PREFIX: Prefix output
'''
########################################################################
#CONSTANT
CWD=os.getcwd()
RSCRIPT=CWD+"/DrawGraphDensity.r"

SUM="Sum"
IDENTIFIED="Identified"
UNIDENTIFIED="Unidentified"
UNASSIGNED="Unassigned"
VIRUSES="Viruses"

UNASSEMBLED="Unassembled"
OTHERS="Others"

STATBYFAMILY="StatByFamily"

TAG_TSV=".tsv"
TAG_OUTPUT1=".Viruses"+TAG_TSV
TAG_OUTPUT2=".Global"+TAG_TSV
EMPTY=""

LIST_SUPP=[UNASSEMBLED,UNIDENTIFIED,OTHERS,VIRUSES]

########################################################################
#Options
parser = OptionParser(conflict_handler="resolve")
parser.add_option("-a","--assembly_stat", dest="assembly_stat")
parser.add_option("-i","--identification_stat", dest="identification_stat")
parser.add_option("-f","--folder_stat", dest="folder_stat")
parser.add_option("-o","--output", dest="output")

(options, args) = parser.parse_args()

sAssembly=options.assembly_stat
if not sAssembly:
	exit("Error : no assembly_stat -a defined, process broken") 

sIdentification=options.identification_stat
if not sIdentification:
	exit("Error : no identification_stat -i defined, process broken")

sFolder=options.folder_stat
if not sFolder:
	exit("Error : no folder_stat -f defined, process broken")
    
sOutput=options.output
if not sOutput:
	exit("Error : no output -o defined, process broken")
    
########################################################################
#Function
def ParseAssemblyStat(sPath):
    dDict={}
    with open(sPath) as oInput:
        bHeader=True
        for sNewLine in oInput:
            if bHeader:
                bHeader=False
                continue
            sLine=sNewLine.strip()
            if len(sLine)==0:
                continue
            tLine=sLine.split("\t")
            sSampleId=tLine[0]
            iUnassembled=int(tLine[-1])
            dDict[sSampleId]=iUnassembled
    return dDict

def ParseIdentificationStat(sPath):
    dDict={}
    with open(sPath) as oInput:
        bHeader=True
        for sNewLine in oInput:
            if bHeader:
                bHeader=False
                continue
            sLine=sNewLine.strip()
            if len(sLine)==0:
                continue
            tLine=sLine.split("\t")
            sSampleId=tLine[0]
            iUnidentified=int(tLine[-1])
            iIdentified=int(tLine[-2])
            dDict[sSampleId]={UNIDENTIFIED:iUnidentified,IDENTIFIED:iIdentified}
    return dDict

def ParseFolderStat(sPath):
    dDict={}
    tSpeciesList=os.listdir(sPath)
    for sSpeciesFile in tSpeciesList:
        sSpeciesName=sSpeciesFile.replace(TAG_TSV,EMPTY)
        dDict[sSpeciesName]={}
        with open(sPath+"/"+sSpeciesFile) as oInput:
            bHeader=True
            for sNewLine in oInput:
                if bHeader:
                    bHeader=False
                    continue
                sLine=sNewLine.strip()
                if len(sLine)==0:
                    continue
                tLine=sLine.split("\t")
                sSampleId=tLine[0]
                iReads=int(tLine[1])
                dDict[sSpeciesName][sSampleId]=iReads
    return dDict

def RegroupVirusData(dVirus):
    dDict={}    
    #Viruses
    for sVirus in dVirus:
        dDict[sVirus]={}
        for sSampleId in dVirus[sVirus]:
            dDict[sVirus][sSampleId]=dVirus[sVirus][sSampleId]
    
    return dDict
    
def WriteRdataframe(dDict,tTarget,sPath):
    with open(sPath,"w") as oOutput:
        sHeader="Legend\tSampleId\tReads\n"
        oOutput.write(sHeader)
        for sKey in tTarget:
            for sSampleId in dDict[sKey]:
                sLine="{}\t{}\t{}\n".format(sKey,sSampleId,dDict[sKey][sSampleId])
                oOutput.write(sLine)

def RegroupGlobalData(dUnassembled,dIdentification,dVirus):
    dVirusReads={}
    for sSampleId in dUnassembled:
        iVirusReads=0
        for sSpecies in dVirus:
            try:
                iVirusReads+=dVirus[sSampleId]
            except KeyError:
                pass
        dVirusReads[sSampleId]=iVirusReads
    
    dDeltaIdentify={}
    for sSampleId in dVirusReads:
        iVirus=dVirusReads[sSampleId]
        iAll=dIdentification[sSampleId][IDENTIFIED]
        iDelta=iAll-iVirus
        dDeltaIdentify[sSampleId]=iDelta
    
    dDict={}    
    #UNASSEMBLED
    dDict[UNASSEMBLED]={}
    for sSampleId in dUnassembled:
        dDict[UNASSEMBLED][sSampleId]=dUnassembled[sSampleId]
    #UNIDENTIFIED
    dDict[UNIDENTIFIED]={}
    for sSampleId in dIdentification:
        dDict[UNIDENTIFIED][sSampleId]=dIdentification[sSampleId][UNIDENTIFIED]
    #OTHERS
    dDict[OTHERS]={}
    for sSampleId in dDeltaIdentify:
        dDict[OTHERS][sSampleId]=dDeltaIdentify[sSampleId]
    #VIRUSES
    dDict[VIRUSES]={}
    for sVirus in dVirus:
        for sSampleId in dVirus[sVirus]:
            try:
                dDict[VIRUSES][sSampleId]+=dVirus[sVirus][sSampleId]
            except KeyError:
                dDict[VIRUSES][sSampleId]=dVirus[sVirus][sSampleId]
    
    return dDict

def LaunchRScript(sPath):
    print("Launch R script...")
    tCommand=["Rscript",RSCRIPT,os.getcwd(),sPath.replace(TAG_TSV,EMPTY)]
    print(" ".join(tCommand)+"$")
    subprocess.call(tCommand)

########################################################################
#MAIN
if __name__ == "__main__":
    dSample2Unassembled=ParseAssemblyStat(sAssembly)
    # print("dSample2Unassembled",dSample2Unassembled)
    dSample2Identification=ParseIdentificationStat(sIdentification)
    # print("dSample2Identification",dSample2Identification)
    dViruses2SampleReads=ParseFolderStat(sFolder)
    # print("dViruses2SampleReads",dViruses2SampleReads)
    dVirusDensity=RegroupVirusData(dViruses2SampleReads)
    dGlobalDensity=RegroupGlobalData(dSample2Unassembled,dSample2Identification,dViruses2SampleReads)
    # print("dDensity",dDensity)
    WriteRdataframe(dVirusDensity,list(sorted(dViruses2SampleReads.keys())),sOutput+TAG_OUTPUT1)
    WriteRdataframe(dGlobalDensity,LIST_SUPP,sOutput+TAG_OUTPUT2)
    LaunchRScript(sOutput+TAG_OUTPUT1)
    LaunchRScript(sOutput+TAG_OUTPUT2)

########################################################################    
iTime2=time.time()
iDeltaTime=iTime2-iTime1
print("Script done: "+str(iDeltaTime))

