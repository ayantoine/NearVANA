# coding: utf-8
"""Python3.6"""
# compatibility: python2.7, python2.6

import time
from optparse import OptionParser

sCurrentVersionScript="v3"
iTime1=time.time()
########################################################################
'''
V3-2019/10/22
Keep IndexEnd value, the index to the end of the linker signals

V2-2019/10/22
No more header in dodeca/midfile
V1-2019/07/01
Create list of specific Kmer from MID_file

python CreateKmerList.py -m MIDFILE -o OUTPUTFILE
MIDFILE: File that contains sample id and associated linker
OUTPUTFILE: Path to the output tab file generated
'''
########################################################################
#CONSTANT

########################################################################
#Options
parser = OptionParser()
parser.add_option("-m","--midfile", dest="midfile")
parser.add_option("-o","--output", dest="output")
parser.add_option("-p","--pid", dest="pid")

(options, args) = parser.parse_args()

sMidFile=options.midfile
if not sMidFile:
	exit("Error : no midfile -m defined, process broken")

sOutput=options.output
if not sOutput:
	exit("Error : no output -o defined, process broken")

sPid=options.pid
if not sPid:
	exit("Error : no pid -p defined, process broken")

########################################################################
#Function 	
def ReadLinkerFile(sString,sPrefix):
	dResult={}
	# bHeader=True
	for sNewLine in open(sString):
		# if bHeader:
			# bHeader=False
			# continue
		sLine=sNewLine.strip()
		tLine=sLine.split("\t")
		print(tLine)
		sRef=sPrefix+"_"+tLine[0]
		sSeq="".join(tLine[1:])
		dResult[sRef]=sSeq
		# print(sRef,sSeq)
	return dResult

def GetKmerRef(dDict):
	iNbrRef=len(dDict)
	iMaxSeqSize=len(dDict[list(dDict.keys())[0]])
	dResult={}
	#For each kmer size avalaible
	for iKmerSize in reversed(range(0,iMaxSeqSize+1)):
		print("Kmer",iKmerSize)
		dResult[iKmerSize]={}
		dKmer2Ref={}
		#If kmer size combination inferior than NbrRef, pass
		if 4**iKmerSize<iNbrRef:
			del dResult[iKmerSize]
			print("Kmer too short to allow at least one kmer-specific by sample, reject library")
			continue
		#For each sequence, store each kmer
		for sSeqId in dDict:
			for i in range(0,iMaxSeqSize+1-iKmerSize):
				iEndIndex=iMaxSeqSize-(i+iKmerSize)
				# print(dDict[sSeqId][i:i+iKmerSize],dDict[sSeqId],iEndIndex)
				try:
					dKmer2Ref[dDict[sSeqId][i:i+iKmerSize]].append((sSeqId,iEndIndex))
				except KeyError:
					dKmer2Ref[dDict[sSeqId][i:i+iKmerSize]]=[(sSeqId,iEndIndex)]
		print("Kmer present",len(dKmer2Ref))
		iKeep=0
		#Keep only kmer that are specific
		for sKmer in dKmer2Ref:
			if len(dKmer2Ref[sKmer])==1:
				dResult[iKmerSize][sKmer]=dKmer2Ref[sKmer][0]
				iKeep+=1
		print("Kmer specific",iKeep)
		#If no kmer specific, break
		if len(dResult[iKmerSize])==0:
			del dResult[iKmerSize]
			print("No specific kmer, reject library")
			break
		#If kmer specific don't cover all ref, remove
		dTemp={}
		for dbValues in dResult[iKmerSize].values():
			dTemp[dbValues[0]]=1
		if len(dTemp)!=iNbrRef:
			del dResult[iKmerSize]
			print("Sample without specific-kmer, reject library")
	
	return dResult

def WriteKmerList(dDict,sPath):
	FILE=open(sPath,"w")
	for iValue in dDict:
		for sString in dDict[iValue]:
			FILE.write(sString+"\t"+dDict[iValue][sString][0]+"\t"+str(dDict[iValue][sString][1])+"\n")
	FILE.close()
	
########################################################################
#MAIN
if __name__ == "__main__":
	dRef2Linker=ReadLinkerFile(sMidFile,sPid)
	dKmerRef=GetKmerRef(dRef2Linker)
	WriteKmerList(dKmerRef,sOutput)
	
########################################################################    
iTime2=time.time()
iDeltaTime=iTime2-iTime1
print("Script done: "+str(iDeltaTime))

