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
UNDERSCORE="_"

ASSEMBLY_BEFORE="Before"
ASSEMBLY_ASSEMBLED="Assembled"
ASSEMBLY_UNASSEMBLED="Unassembled"

########################################################################
#Options
parser = OptionParser()
parser.add_option("-0","--r0", dest="r0")
parser.add_option("-1","--r1", dest="r1")
parser.add_option("-2","--r2", dest="r2")
parser.add_option("-d","--distribution", dest="distribution")
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

sDistribution=options.distribution
if not sDistribution:
    exit("Error : no distribution -d defined, process broken")

sUnmapped=options.unmapped
if not sUnmapped:
    exit("Error : no unmapped -u defined, process broken")

sOutput=options.output
if not sOutput:
    exit("Error : no output -o defined, process broken")
    

########################################################################
#Function     
def GetSample(sPath):
    dDict={}
    for sNewLine in open(sPath):
        sLine=sNewLine.strip()
        if len(sLine)==0:
            continue
        tLine=sNewLine.split(TABULATION)
        sSampleId=tLine[0]
        dDict[sSampleId]={ASSEMBLY_BEFORE:0,ASSEMBLY_ASSEMBLED:0,ASSEMBLY_UNASSEMBLED:0}
    return dDict
        
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
    dSample2Quantity=GetSample(sDistribution)
    dSample2Quantity=CountQuantityBefore(sR0,sR1,sR2,dSample2Quantity)
    dSample2Quantity=CountQuantityRejected(sUnmapped,dSample2Quantity)
    dSample2Quantity=DeductQuantityAssembled(dSample2Quantity)
    WriteOutput(dSample2Quantity,sOutput)
    
########################################################################    
iTime2=time.time()
iDeltaTime=iTime2-iTime1
print("Script done: "+str(iDeltaTime))
