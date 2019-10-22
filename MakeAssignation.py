# coding: utf-8
"""Python3.6"""
# compatibility: python2.7, python2.6

import time
from optparse import OptionParser
import codecs

sCurrentVersionScript="v1"
iTime1=time.time()
########################################################################
'''
V1-2019/07/01
Research Kmer in sequence to assign to a sample

python MakeAssignation.py -1 R1 -2 R2 -k KMERLIST -d WORKDIR -t TAGFILE
R1: Pair-end R1 fastq
R2: Pair-end R2 fastq
KMERLIST: Tabulated file that contains Kmer to sampleId
WORKDIR: Directory where write results
TAGFILE: path of thne tag file, created when correctly ending the script
'''
########################################################################
#CONSTANT
HYPER_SUFFIX="Hyper_Identified.tab"
HYPO1_SUFFIX="Hypo_1_Identified.tab"
HYPO2_SUFFIX="Hypo_2_Identified.tab"
AMBIGUOUS1_SUFFIX="Ambiguous_1.tab"
AMBIGUOUS2_SUFFIX="Ambiguous_2.tab"
UNIDENTIFIED_SUFFIX="Unidentified.tab"

SAMPLENONE="..."
INDEXNONE="."

SEARCHWINDOWS_SIZE=30
########################################################################
#Options
parser = OptionParser()
parser.add_option("-1","--fastq1", dest="fastq1")
parser.add_option("-2","--fastq2", dest="fastq2")
parser.add_option("-k","--kmerlist", dest="kmerlist")
parser.add_option("-d","--workdir", dest="workdir")
parser.add_option("-t","--tagfile", dest="tagfile")

(options, args) = parser.parse_args()

sFastq1=options.fastq1
if not sFastq1:
	exit("Error : no fastq1 -1 defined, process broken")

sFastq2=options.fastq2
if not sFastq2:
	exit("Error : no fastq2 -2 defined, process broken")

sFastq=sFastq1[1:]
if sFastq!=sFastq2[1:]:
	exit("Error : fastq1 and fastq2 must share the same alphabetic designation, process broken")

sKmerList=options.kmerlist
if not sKmerList:
	exit("Error : no kmerlist -k defined, process broken")

sWorkDir=options.workdir
if not sWorkDir:
	exit("Error : no workdir -d defined, process broken")

sTagFile=options.tagfile
if not sTagFile:
	exit("Error : no tagfile -t defined, process broken")

sHyperName=sWorkDir+"/"+sFastq+"_"+HYPER_SUFFIX
sHypo1Name=sWorkDir+"/"+sFastq+"_"+HYPO1_SUFFIX
sHypo2Name=sWorkDir+"/"+sFastq+"_"+HYPO2_SUFFIX
sAmbiguous1Name=sWorkDir+"/"+sFastq+"_"+AMBIGUOUS1_SUFFIX
sAmbiguous2Name=sWorkDir+"/"+sFastq+"_"+AMBIGUOUS2_SUFFIX
sUnidentifiedName=sWorkDir+"/"+sFastq+"_"+UNIDENTIFIED_SUFFIX

########################################################################
#Function 	
def LoadKmerFile(sPath):
	dDict={}
	dOtherDict={}
	for sNewLine in open(sPath):
		sLine=sNewLine.strip()
		tLine=sLine.split()
		sKmer=tLine[0]
		iKmer=len(sKmer)
		sSample=tLine[1]
		iEndIndex=int(tLine[2])
		try:
			dDict[iKmer][sKmer]=sSample
		except KeyError:
			dDict[iKmer]={sKmer:sSample}
		dOtherDict[sKmer]=iEndIndex
	return dDict, dOtherDict
	
def ProcessFastq1(dKmer,dEndIndex,sFastq):
	dResult={}
	iCount=0
	
	for sNewLine in open(sFastq):
		sLine=sNewLine.strip()
		iCount+=1
		if iCount%4==1:
			sSeqId=sLine[1:]
			sSeqCommonId=sSeqId.split(" ")[0]
			sSeqSuffix=sSeqId.split(" ")[1]
		if iCount%4==2:
			sSeq=sLine[:SEARCHWINDOWS_SIZE] #Search linker only on the first nt
			sSample,sEndIndex=AssignSample(sSeq,dKmer,dEndIndex)
			dResult[sSeqCommonId]={"SUFFIX":sSeqSuffix, "SAMPLE":sSample, "END":sEndIndex}	
	return dResult

	
def AssignSample(sString,dKmer,dEndIndex):
	dAssignation={}
	sAssignation=SAMPLENONE
	sEndIndex=INDEXNONE
	#For bigger size of kmer to lower
	for iKmerSize in sorted(dKmer, reverse=True):
		#For all specific kmer
		dPresentKmer={}
		for sKmer in dKmer[iKmerSize]:
			#Check if kmer is present
			if sKmer in sString:
				iKmerStart=sString.find(sKmer)
				iKmerEnd=iKmerStart+len(sKmer)-1
				iEndSignal=iKmerEnd+dEndIndex[sKmer]
				try:
					dPresentKmer[dKmer[iKmerSize][sKmer]].append(iEndSignal)
				except KeyError:
					dPresentKmer[dKmer[iKmerSize][sKmer]]=[iEndSignal]
				#Increase Sample weight
				try:
					dAssignation[dKmer[iKmerSize][sKmer]]+=1
				except KeyError:
					dAssignation[dKmer[iKmerSize][sKmer]]=1
		#If many Sample, clean all weigth 1
		if len(dAssignation)>1:
			if min(dAssignation.values())!=max(dAssignation.values()):
				tClean=[]
				for sSampleId in dAssignation:
					if dAssignation[sSampleId]==1:
						tClean.append(sSampleId)
				for sSampleId in tClean:
					del dAssignation[sSampleId]
		#If kmer of only one sample, assign sequence to this sample
		if len(dAssignation)==1:
			sAssignation=list(dAssignation.keys())[0]
			sEndIndex=str(max(dPresentKmer[dKmer[iKmerSize][sKmer]]))
			break
		#Else, pursue with lower Kmer
	return sAssignation,sEndIndex
	
def ProcessFastq2(dKmer,sFastq,dSeq1):
	HYPERPATH=open(sHyperName,"w")
	HYP01PATH=open(sHypo1Name,"w")
	HYP02PATH=open(sHypo2Name,"w")
	AMBIGUOUS1PATH=open(sAmbiguous1Name,"w")
	AMBIGUOUS2PATH=open(sAmbiguous2Name,"w")
	UNIDENTIFIEDPATH=open(sUnidentifiedName,"w")
	
	iCount=0
	for sNewLine in open(sFastq):
		sLine=sNewLine.strip()
		iCount+=1
		if iCount%4==1:
			sSeqId=sLine[1:]
			sSeqCommonId=sSeqId.split(" ")[0]
			sSeqSuffix=sSeqId.split(" ")[1]
		if iCount%4==2:
			sSeq=sLine[:SEARCHWINDOWS_SIZE] #Search linker only on the first nt
			sSample, sEndIndex=AssignSample(sSeq,dKmer,dEndIndex)
			
			#Same assignation
			if sSample==dSeq1[sSeqCommonId]["SAMPLE"]:
				#Same sample
				if sSample!=SAMPLENONE:
					HYPERPATH.write(sSeqCommonId+" "+dSeq1[sSeqCommonId]["SUFFIX"]+"\t"+dSeq1[sSeqCommonId]["SAMPLE"]+"\t"+dSeq1[sSeqCommonId]["END"]+"\n")
					HYPERPATH.write(sSeqId+"\t"+sSample+"\t"+sEndIndex+"\n")
				else:
					#No sample for both
					UNIDENTIFIEDPATH.write(sSeqCommonId+" "+dSeq1[sSeqCommonId]["SUFFIX"]+"\t"+dSeq1[sSeqCommonId]["SAMPLE"]+"\t"+dSeq1[sSeqCommonId]["END"]+"\n")
					UNIDENTIFIEDPATH.write(sSeqId+"\t"+sSample+"\t"+sEndIndex+"\n")
			else:
				#One assignation and one without
				if sSample==SAMPLENONE or dSeq1[sSeqCommonId]["SAMPLE"]==SAMPLENONE:
					HYP01PATH.write(sSeqCommonId+" "+dSeq1[sSeqCommonId]["SUFFIX"]+"\t"+dSeq1[sSeqCommonId]["SAMPLE"]+"\t"+dSeq1[sSeqCommonId]["END"]+"\n")
					HYP01PATH.write(sSeqId+"\t"+sSample+"\t"+sEndIndex+"\n")
					HYP02PATH.write(sSeqCommonId+" "+dSeq1[sSeqCommonId]["SUFFIX"]+"\t"+SAMPLENONE+"\t"+dSeq1[sSeqCommonId]["END"]+"\n")
					HYP02PATH.write(sSeqId+"\t"+SAMPLENONE+"\t"+sEndIndex+"\n")
				else:
					#Both assigned to sample, but not the same
					AMBIGUOUS1PATH.write(sSeqCommonId+" "+dSeq1[sSeqCommonId]["SUFFIX"]+"\t"+dSeq1[sSeqCommonId]["SAMPLE"]+"\t"+dSeq1[sSeqCommonId]["END"]+"\n")
					AMBIGUOUS1PATH.write(sSeqId+"\t"+sSample+"\t"+sEndIndex+"\n")
					AMBIGUOUS2PATH.write(sSeqCommonId+" "+dSeq1[sSeqCommonId]["SUFFIX"]+"\t"+SAMPLENONE+"\t"+dSeq1[sSeqCommonId]["END"]+"\n")
					AMBIGUOUS2PATH.write(sSeqId+"\t"+SAMPLENONE+"\t"+sEndIndex+"\n")
			
	HYPERPATH.close()
	HYP01PATH.close()
	HYP02PATH.close()
	AMBIGUOUS1PATH.close()
	AMBIGUOUS2PATH.close()
	UNIDENTIFIEDPATH.close()

def CreateTag(sName):
	FILE=open(sName,"w")
	FILE.close()
		
########################################################################
#MAIN
if __name__ == "__main__":
	dKmerRef, dKmerEndIndex=LoadKmerFile(sKmerList)
	dSeqId2Sample=ProcessFastq1(dKmerRef,dKmerEndIndex,sWorkDir+"/"+sFastq1)
	ProcessFastq2(dKmerRef,dKmerEndIndex,sWorkDir+"/"+sFastq2,dSeqId2Sample)
	CreateTag(sTagFile)

########################################################################    
iTime2=time.time()
iDeltaTime=iTime2-iTime1
print("Script done: "+str(iDeltaTime))

