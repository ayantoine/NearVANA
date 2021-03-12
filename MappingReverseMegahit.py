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
V5-2021/03/12
Too much data stored for non-paired data. Try another way.

V4-2020/11/02
Previous version too much long. Try another way, ligther AND faster
V3-2020/10/28
RNAseq data use too much memory. Need a lighter version
V2-2020/10/09
Rework Reads and Pairs comptability
V1-2019/11/04
Extract Sequences from Mapping and dispatch them among file.

python MappingReverseMegahit.py -i INPUT -p PID -m MULTIPLEX
INPUT: path of the mapping file
PID: Project id
MULTIPLEX: Boolean, true if data are multiplexed, false elsewhere
'''
########################################################################
#CONSTANT
REJECT_TAG="@"
NO_MAPPING="*"
COL_SEQID=0 #! count col from 0
COL_MAPPING=2 

CONTIG_BASENAME="Contigs_Megahit"

TRUE=["TRUE","True","true"]
FALSE=["FALSE","False","false"]

UNDERSCORE="_"

SAMPLE="Sample"
READS="Reads"

########################################################################
#Options
parser = OptionParser()
parser.add_option("-i","--input", dest="input")
parser.add_option("-p","--pid", dest="pid")
parser.add_option("-m","--multiplex", dest="multiplex")

(options, args) = parser.parse_args()

sInput=options.input
if not sInput:
	exit("Error : no input -i defined, process broken")

sPID=options.pid
if not sPID:
	exit("Error : no pid -p defined, process broken")

sMultiplex=options.multiplex
if not sMultiplex:
	exit("Error : no multiplex -m defined, process broken")
if sMultiplex in TRUE:
	bMultiplex=True
elif sMultiplex in FALSE:
	bMultiplex=False
else:
	exit("Error : multiplex -m must be one of following: "+",".join(TRUE+FALSE))

#Half-Constant
MEGAHIT_ASSEMBLY_OUTPUT=sPID+"_Temp.Megahit_contigs.fa"

OUTPUT_REVERSEASSEMBLY=sPID+"_All.Megahit_reverseAssembly.tsv"
OUTPUT_MEGAHIT_CONTIG=sPID+"_All.Megahit_contigs.fa"
OUTPUT_R1_UNITIG=sPID+"_R1.Megahit_unassembled.fastq"
OUTPUT_R2_UNITIG=sPID+"_R2.Megahit_unassembled.fastq"
OUTPUT_R0_UNITIG=sPID+"_R0.Megahit_unassembled.fastq"

R0_BASE_FASTQ=sPID+"_R0.Substracted.fastq"
R1_BASE_FASTQ=sPID+"_R1.Substracted.fastq"
R2_BASE_FASTQ=sPID+"_R2.Substracted.fastq"

OUTPUT_REJECTED_CONTIG=sPID+"_All.Megahit_rejectedContigs.fa"
OUTPUT_AMBIGOUS_READ=sPID+"_All.Megahit_ambigousReads.tsv"
OUTPUT_UNMAPPED_READ=sPID+"_All.Megahit_unmappedReads.tsv"

dBASEtoTARGET={
		R0_BASE_FASTQ:OUTPUT_R0_UNITIG,
		R1_BASE_FASTQ:OUTPUT_R1_UNITIG,
		R2_BASE_FASTQ:OUTPUT_R2_UNITIG
	}

print("MEGAHIT_ASSEMBLY_OUTPUT:",MEGAHIT_ASSEMBLY_OUTPUT)
print("OUTPUT_REVERSEASSEMBLY:",OUTPUT_REVERSEASSEMBLY)
print("OUTPUT_MEGAHIT_CONTIG:",OUTPUT_MEGAHIT_CONTIG)
print("OUTPUT_R1_UNITIG:",OUTPUT_R1_UNITIG)
print("OUTPUT_R2_UNITIG:",OUTPUT_R2_UNITIG)
print("OUTPUT_R0_UNITIG:",OUTPUT_R0_UNITIG)
print("R0_BASE_FASTQ:",R0_BASE_FASTQ)
print("R1_BASE_FASTQ:",R1_BASE_FASTQ)
print("R2_BASE_FASTQ:",R2_BASE_FASTQ)
print("OUTPUT_REJECTED_CONTIG:",OUTPUT_REJECTED_CONTIG)
print("OUTPUT_AMBIGOUS_READ:",OUTPUT_AMBIGOUS_READ)
print("OUTPUT_UNMAPPED_READ:",OUTPUT_UNMAPPED_READ)
print("dBASEtoTARGET:",dBASEtoTARGET)

########################################################################
#Function 	
def WriteContigFiles(dContig2Read,bMultiplex):
	FILE_REJECTED=open(OUTPUT_REJECTED_CONTIG,"w")
	FILE_ORIGIN=open(OUTPUT_REVERSEASSEMBLY,"w")
	FILE_FASTA=open(OUTPUT_MEGAHIT_CONTIG,"w")
	
	sId=""
	sShortId=""
	sContent=""
	tSample=[]
	bRejected=False
	iNameId=0
	
	iLineCounter=0
	iTimeStart=time.time()
	
	for sNewLine in open(MEGAHIT_ASSEMBLY_OUTPUT):
		iLineCounter+=1
		if iLineCounter%1000==0:
			iTimeCurrent=time.time()
			iDelta=iTimeCurrent-iTimeStart
			print("Reading line "+str(iLineCounter)+"...\t"+str(iDelta))
			iTimeStart=time.time()
		
		if sNewLine[0]==">":
			if sId!="" and sContent!="":
				if bRejected:
					FILE_REJECTED.write("{}{}".format(sId,sContent))
				else:
					iNameId+=1
					sNewName=CONTIG_BASENAME+"_"+str(iNameId)+"_("+str(len(dContig2Read[sShortId][READS]))+")"
					if bMultiplex:
						tSample=dContig2Read[sShortId][SAMPLE]
					FILE_FASTA.write(">{}\n{}".format(sNewName,sContent))
					for sReadId in dContig2Read[sShortId][READS]:
						FILE_ORIGIN.write(sReadId+"\t"+sNewName+"\t"+sShortId+"\t"+"-".join(sorted(tSample))+"\n")
			sId=sNewLine
			sShortId=sNewLine[1:-1].split(" ")[0]
			sContent=""
			try:
				oCrash=dContig2Read[sShortId]
				bRejected=False
			except KeyError:
				bRejected=True
		else:
			sContent+=sNewLine
	FILE_FASTA.close()
	FILE_ORIGIN.close()
	FILE_REJECTED.close()

def ParseSamfile(sPath):
	# print("Parsing "+sPath)
	dContig2Read={}
	
	FILE_UNMAPPED=open(OUTPUT_UNMAPPED_READ,"w")
	FILE_AMBIGOUS=open(OUTPUT_AMBIGOUS_READ,"w")
	
	iLineCounter=0
	iTimeStart=time.time()

	sPreviousPairId=""
	dReads2Contigs={}
	
	iMappingLine=0
	iUnmapped=0
	iAmbigous=0

	for sNewLine in open(sPath):
		iLineCounter+=1
		if iLineCounter%50000==0:
			iTimeCurrent=time.time()
			iDelta=iTimeCurrent-iTimeStart
			print("Reading line "+str(iLineCounter)+"...\t"+str(iDelta))
			print("\tiMappingLine:{}\tiUnmapped:{}\tiAmbigous:{}".format(iMappingLine,iUnmapped,iAmbigous))
			iTimeStart=time.time()

		if sNewLine[0]==REJECT_TAG:
			continue	

		iMappingLine+=1

		sLine=sNewLine.strip()
		tLine=sLine.split("\t")
		sMap=tLine[COL_MAPPING]
		sId=tLine[COL_SEQID]
		
		if sMap==NO_MAPPING:
			FILE_UNMAPPED.write("{}\n".format(sId))
			iUnmapped+=1
			continue
		
		sCurrentPairId=sId.split(UNDERSCORE)[0]
		if sPreviousPairId=="":
			sPreviousPairId=sCurrentPairId
		if sCurrentPairId==sPreviousPairId:
			try:
				dReads2Contigs[sId].append(sMap)
			except KeyError:
				dReads2Contigs[sId]=[sMap]
		else:
			bAmbigous=False
			for sRead in dReads2Contigs:
				if len(dReads2Contigs[sRead])>1:
					bAmbigous=True
			if bAmbigous:
				iAmbigous+=1
				for sRead in sorted(dReads2Contigs.keys()):
					FILE_AMBIGOUS.write("{}\t{}\n".format(sRead,",".join(dReads2Contigs[sRead])))
			else:
				for sRead in dReads2Contigs:
					sSample=sRead.split(UNDERSCORE)[-1]
					try:
						oCrash=dContig2Read[dReads2Contigs[sRead][0]]
					except KeyError:
						dContig2Read[dReads2Contigs[sRead][0]]={}
					try:
						dContig2Read[dReads2Contigs[sRead][0]][READS].append(sRead)
					except KeyError:
						dContig2Read[dReads2Contigs[sRead][0]][READS]=[sRead]
					try:
						dContig2Read[dReads2Contigs[sRead][0]][SAMPLE][sSample]=True
					except KeyError:
						dContig2Read[dReads2Contigs[sRead][0]][SAMPLE]={sSample:True}
			dReads2Contigs={}
			dReads2Contigs[sId]=[sMap]
		sPreviousPairId=sCurrentPairId
	
	FILE_UNMAPPED.close()
	FILE_AMBIGOUS.close()
	
	return dContig2Read

########################################################################
#MAIN
if __name__ == "__main__":
	dContig2Read=ParseSamfile(sInput)
	# CheckAmbigous(dContig2Read)
	WriteContigFiles(dContig2Read,bMultiplex)

########################################################################    
iTime2=time.time()
iDeltaTime=iTime2-iTime1
print("Script done: "+str(iDeltaTime))
