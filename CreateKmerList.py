# coding: utf-8
"""Python3.6"""
# compatibility: python2.7, python2.6

import time
from optparse import OptionParser

sCurrentVersionScript="v1"
iTime1=time.time()
########################################################################
'''
V1-2019/10/22
No more header in midfile

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
		sRef=sPrefix+"_"+tLine[0]
		sSeq=tLine[1]+tLine[2]
		dResult[sRef]=sSeq
		print(sRef,sSeq)
	return dResult

def GetKmerRef(dDict):
	iNbrRef=len(dDict)
	iMaxSeqSize=len(dDict[list(dDict.keys())[0]])
	dResult={}
	#For each kmer size avalaible
	for iKmerSize in range(0,iMaxSeqSize+1):
		dResult[iKmerSize]={}
		dKmer2Ref={}
		#If kmer size combination inferior than NbrRef, pass
		if 4**iKmerSize<iNbrRef:
			del dResult[iKmerSize]
			continue
		#For each sequence, store each kmer
		for sSeqId in dDict:
			for i in range(0,iMaxSeqSize+1-iKmerSize):
				try:
					dKmer2Ref[dDict[sSeqId][i:i+iKmerSize]].append(sSeqId)
				except KeyError:
					dKmer2Ref[dDict[sSeqId][i:i+iKmerSize]]=[sSeqId]
		#Keep only kmer that are specific
		for sKmer in dKmer2Ref:
			if len(dKmer2Ref[sKmer])==1:
				dResult[iKmerSize][sKmer]=dKmer2Ref[sKmer][0]
		#If no kmer specific, break
		if len(dResult[iKmerSize])==0:
			del dResult[iKmerSize]
			break
		#If kmer specific don't cover all ref, remove
		if len(set(dResult[iKmerSize].values()))!=iNbrRef:
			del dResult[iKmerSize]
	
	return dResult

def WriteKmerList(dDict,sString):
	FILE=open(sString,"w")
	for iValue in dDict:
		for sString in dDict[iValue]:
			FILE.write(sString+"\t"+dDict[iValue][sString]+"\n")
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

