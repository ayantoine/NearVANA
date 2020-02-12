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
########################################################################
#Options
parser = OptionParser()
parser.add_option("-f","--fastq", dest="fastq")
parser.add_option("-r","--ref", dest="ref")
parser.add_option("-s","--sample", dest="sample")
parser.add_option("-i","--pairid", dest="pairid")

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

########################################################################
#Function 	
def LoadRef(sPath,sSample,sPair):
	dResult={}
	for sNewLine in open(sPath):
		if sSample in sNewLine:
			if sPair+ILLUMINA_PAIR_TAG in sNewLine:
				sLine=sNewLine.strip()
				tLine=sLine.split("\t")
				sSeqName=tLine[0]
				iEndIndex=tLine[2]
				dResult[sSeqName]=int(iEndIndex)
	return dResult

def WriteSplitFastq(sPath,dList,sSID):
	FILE=open(sSID+"/"+sSID+"_"+sPath+"."+SPLIT_TAG,"w")
	sSeqName=""
	sContent=""
	sInterline=""
	sQuality=""
	iLineCounter=0
	iSeqAssociated=0
	for sNewLine in open(sPath):
		iLineCounter+=1
		if iLineCounter%4==1:
			if sSeqName!="":
				try:
					iEndIndex=dList[sSeqName[1:-1]] #remove starting @ and ending \n
					# print(iEndIndex)
					# print(sSeqName)
					# print(sContent)
					# print(sInterline)
					# print(sQuality)
					sSeqName=sSeqName.replace("\n"," "+sSID+"\n") 
					FILE.write(sSeqName+sContent[iEndIndex:]+sInterline+sQuality[iEndIndex:])
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
			FILE.write(sSeqName+sContent[iEndIndex:]+sInterline+sQuality[iEndIndex:])
			iSeqAssociated+=1
		except KeyError:
			pass
	
	print(sSID+"/"+sSID+"_"+sPath+"."+SPLIT_TAG+" contains "+str(iSeqAssociated)+" sequences")
	FILE.close()
		
########################################################################
#MAIN
if __name__ == "__main__":
	dListOfSeq=LoadRef(sRef,sSampleTag,sPairId)
	WriteSplitFastq(sFastq,dListOfSeq,sSampleId)
	
	
########################################################################    
iTime2=time.time()
iDeltaTime=iTime2-iTime1
print("Script done: "+str(iDeltaTime))

