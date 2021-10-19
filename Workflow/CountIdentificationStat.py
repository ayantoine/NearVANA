# coding: utf-8
"""Python3.6"""
# compatibility: python2.7, python2.6

import time
from optparse import OptionParser

sCurrentVersionScript="v1"
iTime1=time.time()
########################################################################
'''
V1-2021/09/24
Count reads repartition by sample before asembly, after assembly and unassembled

python CountIdentificationStat.py -t TABLE -s ASSEMBLY_STAT -o OUTPUT
TABLE: NearVANA results tab
ASSEMBLY_STAT: CountAssemblyStat.py output
OUTPUT: Final file
'''
########################################################################
#CONSTANT
DEFAULT="UnassignedReads"

TABULATION="\t"
UNDERSCORE="_"

IDENTIFICATION_BEFORE="Before"
IDENTIFICATION_IDENTIFIED="Identified"
IDENTIFICATION_UNIDENTIFIED="NotIdentified"

########################################################################
#Options
parser = OptionParser()
parser.add_option("-t","--table", dest="table")
parser.add_option("-s","--stat", dest="stat")
parser.add_option("-o","--output", dest="output")

(options, args) = parser.parse_args()

sTable=options.table
if not sTable:
    exit("Error : no table -t defined, process broken")

sStat=options.stat
if not sStat:
    exit("Error : no stat -s defined, process broken")

sOutput=options.output
if not sOutput:
    exit("Error : no output -o defined, process broken")
    

########################################################################
#Function     
def GetSampleFromStat(sPath):
    dDict={}
    bHeader=True
    for sNewLine in open(sPath):
        if bHeader:
            bHeader=False
            continue
        sLine=sNewLine.strip()
        if len(sLine)==0:
            continue
        tLine=sNewLine.split(TABULATION)
        sSampleId=tLine[0]
        iValue=int(tLine[2])
        dDict[sSampleId]={IDENTIFICATION_BEFORE:iValue,IDENTIFICATION_IDENTIFIED:0,IDENTIFICATION_UNIDENTIFIED:0}
    return dDict
        
def CountQuantityIdentified(sPath,dDict):
    bHeader=True
    for sNewLine in open(sPath):
        if bHeader:
            bHeader=False
            continue
        sLine=sNewLine.strip()
        if len(sLine)==0:
            continue
        tLine=sNewLine.split(TABULATION)
        sSampleId=tLine[2]
        iValue=int(tLine[3])
        dDict[sSampleId][IDENTIFICATION_IDENTIFIED]+=iValue
    return dDict

def DeductQuantityRejected(dDict):
    for sKey in dDict:
        dDict[sKey][IDENTIFICATION_UNIDENTIFIED]=dDict[sKey][IDENTIFICATION_BEFORE]-dDict[sKey][IDENTIFICATION_IDENTIFIED]
    return dDict

def WriteOutput(dDict,sOutput):
    FILE=open(sOutput,"w")
    sHeader="{}\t{}\t{}\t{}\n".format("SampleId","Reads before","Reads identified","Reads unidentified")
    FILE.write(sHeader)
    for sKey in sorted(dDict):
        sLine="{}\t{}\t{}\t{}\n".format(sKey,dDict[sKey][IDENTIFICATION_BEFORE],
                                        dDict[sKey][IDENTIFICATION_IDENTIFIED],
                                        dDict[sKey][IDENTIFICATION_UNIDENTIFIED])
        FILE.write(sLine)
    FILE.close()

########################################################################
#MAIN
if __name__ == "__main__":
    dSample2Quantity=GetSampleFromStat(sStat)
    dSample2Quantity=CountQuantityIdentified(sTable,dSample2Quantity)
    dSample2Quantity=DeductQuantityRejected(dSample2Quantity)
    WriteOutput(dSample2Quantity,sOutput)
    
########################################################################    
iTime2=time.time()
iDeltaTime=iTime2-iTime1
print("Script done: "+str(iDeltaTime))
