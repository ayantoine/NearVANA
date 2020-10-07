# coding: utf-8
"""Python3.6"""
# compatibility: python2.7, python2.6

import time
from optparse import OptionParser

sCurrentVersionScript="v3"
iTime1=time.time()
########################################################################
'''
V3-2020/02/12
Adapt to multiplate anaysis

V2-2019/10/23
Add SampleId into name
V1-2019/10/22
Split fastq by SampleId and store into specific SampleFolder (All SampleFolder
must be already existing)
Remove linker during process

python SplitReads.py -f FASTQ -r REFFILE -s SAMPLE -p PID -i PAIRID
FASTQ: Fastq with all sequences
REFFILE: ${PID}_Hyper_Identified.tsv
SAMPLE: SampleId
PID: Processus Id
PAIRID: id of the pair end file
'''
########################################################################
#CONSTANT
ILLUMINA_PAIR_TAG=":N:0:"

SPLIT_TAG="split"

UNASSIGNED="UnassignedReads"
NO_SAMPLE="..."
NO_INDEX="."

DEBUG=True
########################################################################
#Options
parser = OptionParser()
parser.add_option("-f","--fastq", dest="fastq")
parser.add_option("-r","--ref", dest="ref")
parser.add_option("-s","--sample", dest="sample")
parser.add_option("-i","--pairid", dest="pairid")
parser.add_option("-o","--output", dest="output")
parser.add_option("-u","--use_unassigned", dest="use_unassigned")

(options, args) = parser.parse_args()

sFastq=options.fastq
if not sFastq:
	exit("Error : no fastq -f defined, process broken")

sRef=options.ref
if not sRef:
	exit("Error : no ref -r defined, process broken")
	
sPairId=options.pairid
if not sPairId:
	exit("Error : no pairid -i defined, process broken")

sSampleId=options.sample
if not sSampleId:
	exit("Error : no sample -s defined, process broken")

sOutput=options.output
if not sOutput:
	exit("Error : no output -o defined, process broken")

sUnassigned=options.use_unassigned
if not sUnassigned:
	bUnassigned=False
else:
	bUnassigned=True

########################################################################
#Function 	
def LoadRef(sPath,sSample,sPair,bKeepUnassigned):
	dResult={}
	if bKeepUnassigned and sSample==UNASSIGNED:
		sSample=NO_SAMPLE
	for sNewLine in open(sPath):
		if sSample in sNewLine:
			if sPair+ILLUMINA_PAIR_TAG in sNewLine:
				sLine=sNewLine.strip()
				tLine=sLine.split("\t")
				sSeqName=tLine[0]
				sEndIndex=tLine[2]
				if sEndIndex==NO_INDEX:
					sEndIndex="0"
				dResult[sSeqName]=int(sEndIndex)
	return dResult

def WriteSplitFastq(sPath,dList,sSID,sOut):
	if 
	FILE=open(sSID+"/"+sOut,"a")
	sSeqName=""
	sContent=""
	sInterline=""
	sQuality=""
	iLineCounter=0
	iSeqAssociated=0
	iEmptySeq=0
	for sNewLine in open(sPath):
		iLineCounter+=1
		if iLineCounter%4==1:
			if sSeqName!="":
				try:
					iEndIndex=dList[sSeqName[1:-1]] 
					sSeqName=sSeqName.replace("\n"," "+sSID+"\n") 
					sLine=sSeqName+sContent[iEndIndex:]+sInterline+sQuality[iEndIndex:]
					if sLine.count("\n")!=4:
						iEmptySeq+=1
						continue
					FILE.write(sLine)
					iSeqAssociated+=1
				except KeyError:
					pass
			sSeqName=sNewLine
			sContent=""
			sInterline=""
			sQuality=""
		elif iLineCounter%4==2:
			sContent+=sNewLine
		elif iLineCounter%4==3:
			sInterline+=sNewLine
		else:
			sQuality+=sNewLine
	if sSeqName!="":
		try:
			iEndIndex=dList[sSeqName[1:-1]] #remove starting @ and ending \n
			sSeqName=sSeqName.replace("\n"," "+sSID+"\n") 
			sLine=sSeqName+sContent[iEndIndex:]+sInterline+sQuality[iEndIndex:]
			if sLine.count("\n")!=4:
				iEmptySeq+=1
				pass
			else:
				FILE.write(sLine)
				iSeqAssociated+=1
		except KeyError:
			pass
	
	print(sSID+"/"+sSID+"_"+sPath+"."+SPLIT_TAG+" contains "+str(iSeqAssociated)+" sequences")
	print(str(iEmptySeq)+" empty sequences were rejected")
	FILE.close()
		
########################################################################
#MAIN
if __name__ == "__main__":
	dListOfSeq=LoadRef(sRef,sSampleId,sPairId,bUnassigned)
	WriteSplitFastq(sFastq,dListOfSeq,sSampleId,sOutput)
	
	
########################################################################    
iTime2=time.time()
iDeltaTime=iTime2-iTime1
print("Script done: "+str(iDeltaTime))

