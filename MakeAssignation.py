# coding: utf-8
"""Python3.6"""
# compatibility: python2.7, python2.6

import time
from optparse import OptionParser
import codecs

sCurrentVersionScript="v3"
iTime1=time.time()
########################################################################
'''
V3-2020/01/21
Adapt to array wioth limited value

V2-2019/10/30
Work with index on base file instead multiple subfile (decrease memory usage)
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
HYPER_SUFFIX="Hyper_Identified.tsv"
HYPO1_SUFFIX="Hypo_1_Identified.tsv"
HYPO2_SUFFIX="Hypo_2_Identified.tsv"
AMBIGUOUS1_SUFFIX="Ambiguous_1.tsv"
AMBIGUOUS2_SUFFIX="Ambiguous_2.tsv"
UNIDENTIFIED_SUFFIX="Unidentified.tsv"

SAMPLENONE="..."
INDEXNONE="."

DEFAULT_SEQ_BY_TASK=1000000 #Must be similar in QsubAssignation.py
LINE_BY_FASTQ=4
# LINE_BY_TASK=SEQ_BY_TASK*LINE_BY_FASTQ

SEARCHWINDOWS_SIZE=30

CONF_COMMENT="#"
CONF_STEP="="
KEYCONF_SCALL="SCALL"
KEYCONF_SPARAM="SPARAM"
KEYCONF_STASKARRAY="STASKARRAY"
KEYCONF_SMAXTASK="SMAXTASK"
KEYCONF_SMAXSIMJOB="SMAXSIMJOB"
KEYCONF_SMAXARRAYSIZE="SMAXARRAYSIZE"
KEYCONF_STASKID="STASKID"
KEYCONF_SPSEUDOTASKID="SPSEUDOTASKID"
########################################################################
#Options
parser = OptionParser()
parser.add_option("-1","--fastq1", dest="fastq1")
parser.add_option("-2","--fastq2", dest="fastq2")
parser.add_option("-k","--kmerlist", dest="kmerlist")
parser.add_option("-d","--workdir", dest="workdir")
parser.add_option("-t","--tagfile", dest="tagfile")
parser.add_option("-i","--index", dest="index")
parser.add_option("-q","--quantity", dest="quantity")
parser.add_option("-c","--conffile", dest="conffile")

(options, args) = parser.parse_args()

sFastq1=options.fastq1
if not sFastq1:
	exit("Error : no fastq1 -1 defined, process broken")

sFastq2=options.fastq2
if not sFastq2:
	exit("Error : no fastq2 -2 defined, process broken")

sKmerList=options.kmerlist
if not sKmerList:
	exit("Error : no kmerlist -k defined, process broken")

sWorkDir=options.workdir
if not sWorkDir:
	exit("Error : no workdir -d defined, process broken")

sTagFile=options.tagfile
if not sTagFile:
	exit("Error : no tagfile -t defined, process broken")
	
sIndex=options.index
if not sIndex:
	exit("Error : no index -i defined, process broken")
try:
	iIndex=int(sIndex)
except ValueError:
	exit("Error : index -i must be an integer, process broken")
	
sQuantity=options.quantity
if not sQuantity:
	exit("Error : no quantity -q defined, process broken")
try:
	sQuantity=sQuantity.split(".")[0]
	iQuantity=int(sQuantity)
except ValueError:
	exit("Error : quantity -q must be an integer, process broken")

sConf=options.conffile
if not sConf:
	exit("Error : no conffile -c defined, process broken")

sHyperName=sWorkDir+"/"+str(iIndex)+"_"+HYPER_SUFFIX
sHypo1Name=sWorkDir+"/"+str(iIndex)+"_"+HYPO1_SUFFIX
sHypo2Name=sWorkDir+"/"+str(iIndex)+"_"+HYPO2_SUFFIX
sAmbiguous1Name=sWorkDir+"/"+str(iIndex)+"_"+AMBIGUOUS1_SUFFIX
sAmbiguous2Name=sWorkDir+"/"+str(iIndex)+"_"+AMBIGUOUS2_SUFFIX
sUnidentifiedName=sWorkDir+"/"+str(iIndex)+"_"+UNIDENTIFIED_SUFFIX

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
	
def ProcessFastq1(dKmer,dEndIndex,sFastq,iIndex,iJobByTask):
	dResult={}
	iCount=0
	iLineByTask=iJobByTask*LINE_BY_FASTQ
	
	for sNewLine in open(sFastq):
		iCount+=1
		if iCount<iIndex*iLineByTask+1:
			continue
		elif iCount>=(iIndex+1)*iLineByTask+1:
			break
		sLine=sNewLine.strip()
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
	dPresentKmer={}
	sAssignation=SAMPLENONE
	sEndIndex=INDEXNONE
	#For bigger size of kmer to lower
	for iKmerSize in sorted(dKmer, reverse=True):
		#For all specific kmer
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
		# ##DEBUG
		# if len(dAssignation)>2:
			# print(dAssignation)
			# print(dPresentKmer)
		#If many Sample, clean all weigth 1
		if len(dAssignation)>1:
			if min(dAssignation.values())!=max(dAssignation.values()):
				tClean=[]
				for sSampleId in dAssignation:
					if dAssignation[sSampleId]==1:
						tClean.append(sSampleId)
				for sSampleId in tClean:
					del dAssignation[sSampleId]
					del dPresentKmer[sSampleId]
		#If kmer of only one sample, assign sequence to this sample
		if len(dAssignation)==1:
			sAssignation=list(dAssignation.keys())[0]
			sEndIndex=str(max(dPresentKmer[sAssignation]))
			break
		#Else, pursue with lower Kmer
	return sAssignation,sEndIndex
	
def ProcessFastq2(dKmer,dEndIndex,sFastq,dSeq1,iIndex,iJobByTask):
	HYPERPATH=open(sHyperName,"w")
	HYP01PATH=open(sHypo1Name,"w")
	HYP02PATH=open(sHypo2Name,"w")
	AMBIGUOUS1PATH=open(sAmbiguous1Name,"w")
	AMBIGUOUS2PATH=open(sAmbiguous2Name,"w")
	UNIDENTIFIEDPATH=open(sUnidentifiedName,"w")
	
	iCount=0
	iLineByTask=iJobByTask*LINE_BY_FASTQ
	for sNewLine in open(sFastq):
		iCount+=1
		if iCount<iIndex*iLineByTask+1:
			continue
		elif iCount>=(iIndex+1)*iLineByTask+1:
			break
		sLine=sNewLine.strip()
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

def LoadConfFile(sPath):
	dDict={}
	for sNewLine in open(sPath):
		# sLine=sNewLine.strip()
		sLine=sNewLine.replace("\n","")
		sLine=sLine.replace("\\","")
		sLine=sLine.replace('"','')
		if len(sLine)==0:
			continue
		if sLine[0]==CONF_COMMENT:
			continue
		tLine=sLine.split(CONF_STEP)
		dDict[tLine[0]]=CONF_STEP.join(tLine[1:])
	return dDict	

def GetJobByTask(iSeq,iTask):
	if iTask==0:
		return DEFAULT_SEQ_BY_TASK
	if iSeq%iTask==0:
		iResult=iSeq/iTask
	else:
		iResult=int(round(iSeq/iTask,0))+1
	return iResult
		
########################################################################
#MAIN
if __name__ == "__main__":
	dConf=LoadConfFile(sConf)
	iJobByTask=GetJobByTask(iQuantity,int(dConf[KEYCONF_SMAXARRAYSIZE]))
	dKmerRef, dKmerEndIndex=LoadKmerFile(sKmerList)
	dSeqId2Sample=ProcessFastq1(dKmerRef,dKmerEndIndex,sFastq1,iIndex,iJobByTask)
	ProcessFastq2(dKmerRef,dKmerEndIndex,sFastq2,dSeqId2Sample,iIndex,iJobByTask)
	CreateTag(sTagFile)

########################################################################    
iTime2=time.time()
iDeltaTime=iTime2-iTime1
print("Script done: "+str(iDeltaTime))

