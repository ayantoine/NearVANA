# coding: utf-8
"""Python3.6"""
# compatibility: python2.7, python2.6

import time
from optparse import OptionParser
import re

sCurrentVersionScript="v2"
iTime1=time.time()
########################################################################
'''
V2-2021/09/27
Switch on data to get back Sample Id. Distribution can't work for multiplate analysis.

V1-2021/09/24
Count reads repartition by sample before asembly, after assembly and unassembled

python CountAssemblyStat.py -0 R0 -1 R1 -2 R2 -u UNMAPPED -o OUTPUT
DISTRIBUTION: NearVANA Hyper_Distribution.tsv
R0: R0 fastq input
R1: R1 fastq input
R2: R2 fastq input
UNMAPPED: NearVANA Megahit_unmappedReads.tsv
OUTPUT: Final file
'''
########################################################################
#CONSTANT
DEFAULT="UnassignedReads"

TABULATION="\t"
SPACE=" "
UNDERSCORE="_"
EQUAL="="
EMPTY=""

ASSEMBLY_BEFORE="Before"
ASSEMBLY_ASSEMBLED="Assembled"
ASSEMBLY_UNASSEMBLED="Unassembled"

DATA_COMMENT_TAG="#"
DATA_PLATE_TAG="PLATE"
DATA_PARENTHESIS_REGEX="\(|\)"
DATA_TABULATION="\t"

########################################################################
#Options
parser = OptionParser()
parser.add_option("-0","--r0", dest="r0")
parser.add_option("-1","--r1", dest="r1")
parser.add_option("-2","--r2", dest="r2")
parser.add_option("-d","--datafile", dest="datafile")
parser.add_option("-u","--unmapped", dest="unmapped")
parser.add_option("-o","--output", dest="output")

(options, args) = parser.parse_args()

sR0=options.r0
if not sR0:
    exit("Error : no r0 -0 defined, process broken")

sR1=options.r1
if not sR1:
    exit("Error : no r1 -1 defined, process broken")
    
sR2=options.r2
if not sR2:
    exit("Error : no r2 -2 defined, process broken")

sDatafile=options.datafile
if not sDatafile:
    exit("Error : no datafile -d defined, process broken")

sUnmapped=options.unmapped
if not sUnmapped:
    exit("Error : no unmapped -u defined, process broken")

sOutput=options.output
if not sOutput:
    exit("Error : no output -o defined, process broken")
    

########################################################################
#Function     
def Convert2SampleDict(tList):
    dDict={}
    for sKey in tList:
    dDict[sKey]={ASSEMBLY_BEFORE:0,ASSEMBLY_ASSEMBLED:0,ASSEMBLY_UNASSEMBLED:0}
    return dDict

def GetSampleList(sPath):
    dTag2File={}
    tTag=[]
    for sNewLine in open(sPath):
        sLine=sNewLine.strip()
        if len(sLine)==0:
            continue
        if sLine[0]==DATA_COMMENT_TAG:
            continue
        if DATA_PLATE_TAG in sLine:
            sPart2=sLine.split(EQUAL)[-1]
            sContent=EMPTY.join(re.split(DATA_PARENTHESIS_REGEX,sPart2))
            tTag=sContent.split(SPACE)
            continue
        for sTag in tTag:
            if sTag in sLine:
                sPart2=sLine.split(EQUAL)[-1]
                sContent=EMPTY.join(re.split(DATA_PARENTHESIS_REGEX,sLine))
                tContent=sContent.split(SPACE)
                sFile=tContent[-1]
                dTag2File[sTag]=sFile
    tResult=[]
    for sTag in sorted(dTag2File):
        for sNewLine in open(dTag2File[sTag]):
            sLine=sNewLine.strip()
            if len(sLine)==0:
                continue
            tLine=sLine.split(DATA_TABULATION)
            sSampleId=sTag+tLine[0]
            tResult.append(sSampleId)
    return tResult
    
def CountQuantityBefore(sR0,sR1,sR2,dDict):
    for sFastq in [sR0,sR1,sR2]:
        iLineCounter=0
        for sNewLine in open(sFastq):
            iLineCounter+=1
            if iLineCounter%4==1:
                sLine=sNewLine.strip()
                sSampleId=sLine.split(UNDERSCORE)[-1]
                dDict[sSampleId][ASSEMBLY_BEFORE]+=1
    return dDict

def CountQuantityRejected(sUnmapped,dDict):
    for sNewLine in open(sUnmapped):
        sLine=sNewLine.strip()
        sSampleId=sLine.split(UNDERSCORE)[-1]
        dDict[sSampleId][ASSEMBLY_UNASSEMBLED]+=1
    return dDict

def DeductQuantityAssembled(dDict):
    for sKey in dDict:
        dDict[sKey][ASSEMBLY_ASSEMBLED]=dDict[sKey][ASSEMBLY_BEFORE]-dDict[sKey][ASSEMBLY_UNASSEMBLED]
    return dDict

def WriteOutput(dDict,sOutput):
    FILE=open(sOutput,"w")
    sHeader="{}\t{}\t{}\t{}\n".format("SampleId","Reads before","Reads assembled","Reads unassembled")
    FILE.write(sHeader)
    for sKey in sorted(dDict):
        sLine="{}\t{}\t{}\t{}\n".format(sKey,dDict[sKey][ASSEMBLY_BEFORE],
                                        dDict[sKey][ASSEMBLY_ASSEMBLED],
                                        dDict[sKey][ASSEMBLY_UNASSEMBLED])
        FILE.write(sLine)
    FILE.close()

########################################################################
#MAIN
if __name__ == "__main__":
    print("Retrieve all sample Id...")
    tSampleList=GetSampleList(sDatafile)
    dSample2Quantity=Convert2SampleDict(tSampleList)
    dSample2Quantity=CountQuantityBefore(sR0,sR1,sR2,dSample2Quantity)
    dSample2Quantity=CountQuantityRejected(sUnmapped,dSample2Quantity)
    dSample2Quantity=DeductQuantityAssembled(dSample2Quantity)
    WriteOutput(dSample2Quantity,sOutput)
    
########################################################################    
iTime2=time.time()
iDeltaTime=iTime2-iTime1
print("Script done: "+str(iDeltaTime))
