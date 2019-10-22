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
	for sNewLine in open(sPath):
		sLine=sNewLine.strip()
		tLine=sLine.split()
		sKmer=tLine[0]
		iKmer=len(sKmer)
		sSample=tLine[1]
		try:
			dDict[iKmer][sKmer]=sSample
		except KeyError:
			dDict[iKmer]={sKmer:sSample}
	return dDict
	
def ProcessFastq1(dKmer,sFastq):
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
			sSample=AssignSample(sSeq,dKmer)
			dResult[sSeqCommonId]={"SUFFIX":sSeqSuffix, "SAMPLE":sSample}	
	return dResult

	
def AssignSample(sString,dKmer):
	dAssignation={}
	sAssignation=SAMPLENONE
	#For bigger size of kmer to lower
	for iKmerSize in sorted(dKmer, reverse=True):
		#For all specific kmer
		for sKmer in dKmer[iKmerSize]:
			#Check if kmer is present
			if sKmer in sString:
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
			break
		#Else, pursue with lower Kmer
	return sAssignation
	
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
			sSample=AssignSample(sSeq,dKmer)
			
			#Same assignation
			if sSample==dSeq1[sSeqCommonId]["SAMPLE"]:
				#Same sample
				if sSample!=SAMPLENONE:
					HYPERPATH.write(sSeqCommonId+" "+dSeq1[sSeqCommonId]["SUFFIX"]+"\t"+dSeq1[sSeqCommonId]["SAMPLE"]+"\n")
					HYPERPATH.write(sSeqId+"\t"+sSample+"\n")
				else:
					#No sample for both
					UNIDENTIFIEDPATH.write(sSeqCommonId+" "+dSeq1[sSeqCommonId]["SUFFIX"]+"\t"+dSeq1[sSeqCommonId]["SAMPLE"]+"\n")
					UNIDENTIFIEDPATH.write(sSeqId+"\t"+sSample+"\n")
			else:
				#One assignation and one without
				if sSample==SAMPLENONE or dSeq1[sSeqCommonId]["SAMPLE"]==SAMPLENONE:
					HYP01PATH.write(sSeqCommonId+" "+dSeq1[sSeqCommonId]["SUFFIX"]+"\t"+dSeq1[sSeqCommonId]["SAMPLE"]+"\n")
					HYP01PATH.write(sSeqId+"\t"+sSample+"\n")
					HYP02PATH.write(sSeqCommonId+" "+dSeq1[sSeqCommonId]["SUFFIX"]+"\t"+SAMPLENONE+"\n")
					HYP02PATH.write(sSeqId+"\t"+SAMPLENONE+"\n")
				else:
					#Both assigned to sample, but not the same
					AMBIGUOUS1PATH.write(sSeqCommonId+" "+dSeq1[sSeqCommonId]["SUFFIX"]+"\t"+dSeq1[sSeqCommonId]["SAMPLE"]+"\n")
					AMBIGUOUS1PATH.write(sSeqId+"\t"+sSample+"\n")
					AMBIGUOUS2PATH.write(sSeqCommonId+" "+dSeq1[sSeqCommonId]["SUFFIX"]+"\t"+SAMPLENONE+"\n")
					AMBIGUOUS2PATH.write(sSeqId+"\t"+SAMPLENONE+"\n")
			
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
	dKmerRef=LoadKmerFile(sKmerList)
	dSeqId2Sample=ProcessFastq1(dKmerRef,sWorkDir+"/"+sFastq1)
	ProcessFastq2(dKmerRef,sWorkDir+"/"+sFastq2,dSeqId2Sample)
	CreateTag(sTagFile)

########################################################################    
iTime2=time.time()
iDeltaTime=iTime2-iTime1
print("Script done: "+str(iDeltaTime))

