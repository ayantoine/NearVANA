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
V2-2020/10/09
Rework Reads and Pairs comptability

V1-2019/11/04
Extract Sequences from Mapping and dispatch them among file.

python MappingReverseMegahit.py -i INPUT -p PID
INPUT: path of the mapping file
PID: Project id
'''
########################################################################
#CONSTANT
REJECT_TAG="@"
NO_MAPPING="*"
COL_SEQID=0 #! count col from 0
COL_MAPPING=2 

CONTIG_BASENAME="Contigs_Megahit"


########################################################################
#Options
parser = OptionParser()
parser.add_option("-i","--input", dest="input")
parser.add_option("-p","--pid", dest="pid")

(options, args) = parser.parse_args()

sInput=options.input
if not sInput:
	exit("Error : no input -i defined, process broken")

sPID=options.pid
if not sPID:
	exit("Error : no pid -p defined, process broken")

#Half-Constant
MEGAHIT_ASSEMBLY_OUTPUT=sPID+"_Temp.Megahit_contigs.fa"

OUTPUT_FILENAME=sPID+"_All.Megahit_reverseAssembly.tsv"
OUTPUT_SPADES_CONTIG=sPID+"_All.Megahit_contigs.fa"
OUTPUT_R1_UNITIG=sPID+"_R1.Megahit_unassembled.fastq"
OUTPUT_R2_UNITIG=sPID+"_R2.Megahit_unassembled.fastq"
OUTPUT_R0_UNITIG=sPID+"_R0.Megahit_unassembled.fastq"

R0_BASE_FASTQ=sPID+"_R0.Substracted.fastq"
R1_BASE_FASTQ=sPID+"_R1.Substracted.fastq"
R2_BASE_FASTQ=sPID+"_R2.Substracted.fastq"

OUTPUT_REJECTED_CONTIG=sPID+"_All.Megahit_rejectedContigs.fa"
OUTPUT_AMBIGOUS_READ=sPID+"_All.Megahit_ambigousReads.tsv"

dBASEtoTARGET={
		R0_BASE_FASTQ:OUTPUT_R0_UNITIG,
		R1_BASE_FASTQ:OUTPUT_R1_UNITIG,
		R2_BASE_FASTQ:OUTPUT_R2_UNITIG
	}

print("MEGAHIT_ASSEMBLY_OUTPUT:",MEGAHIT_ASSEMBLY_OUTPUT)
print("OUTPUT_FILENAME:",OUTPUT_FILENAME)
print("OUTPUT_SPADES_CONTIG:",OUTPUT_SPADES_CONTIG)
print("OUTPUT_R1_UNITIG:",OUTPUT_R1_UNITIG)
print("OUTPUT_R2_UNITIG:",OUTPUT_R2_UNITIG)
print("OUTPUT_R0_UNITIG:",OUTPUT_R0_UNITIG)
print("R0_BASE_FASTQ:",R0_BASE_FASTQ)
print("R1_BASE_FASTQ:",R1_BASE_FASTQ)
print("R2_BASE_FASTQ:",R2_BASE_FASTQ)
print("OUTPUT_REJECTED_CONTIG:",OUTPUT_REJECTED_CONTIG)
print("OUTPUT_AMBIGOUS_READ:",OUTPUT_AMBIGOUS_READ)
print("dBASEtoTARGET:",dBASEtoTARGET)

########################################################################
#Function 	
# def Parse(sPath):
	# print("Parsing "+sPath)
	# dContig2Read={}
	# dContig2Sample={}
	
	# dPair2Contig={}
	# dPair2Read2Contig={}

	# iLineCounter=0

	# for sNewLine in open(sPath):
		# iLineCounter+=1
		# if iLineCounter%50000==0:
			# print("Reading line "+str(iLineCounter)+"...")
		# if sNewLine[0]==REJECT_TAG:
			# continue
		# sLine=sNewLine.strip()
		# tLine=sLine.split("\t")
		# sMap=tLine[COL_MAPPING]
		# sId=tLine[COL_SEQID]
		
		# try:
			# dContig2Read[sMap].add(sId)
		# except KeyError:
			# dContig2Read[sMap]=set([sId])
			
		# try:
			# dContig2Sample[sMap].add(sId.split("_")[-1])
		# except KeyError:
			# dContig2Sample[sMap]=set([sId.split("_")[-1]])
			
		# try:
			# dPair2Contig[sId.split("_")[0]].add(sMap)
		# except KeyError:
			# dPair2Contig[sId.split("_")[0]]=set([sMap])
		
		# try:
			# dPair2Read2Contig[sId.split("_")[0]][sId]=sMap
		# except KeyError:
			# dPair2Read2Contig[sId.split("_")[0]]={sId:sMap}
	
	# setLost=set([])
	# FILE=open(OUTPUT_AMBIGOUS_READ,"w")
	# iPairOk=0
	# iPairLost=0
	# iPairHalfLost=0
	# iPairAmbigous=0
	# for sKey in dPair2Contig:
		# if NO_MAPPING in dPair2Contig[sKey] and len(dPair2Contig[sKey])==1:
			# iPairLost+=1
			# setLost.add(sKey)
		# elif NO_MAPPING in dPair2Contig[sKey] and len(dPair2Contig[sKey])!=1:
			# iPairHalfLost+=1
		# elif NO_MAPPING not in dPair2Contig[sKey] and len(dPair2Contig[sKey])==1:
			# iPairOk+=1
		# elif NO_MAPPING not in dPair2Contig[sKey] and len(dPair2Contig[sKey])!=1:
			# iPairAmbigous+=1
			# for sRead in dPair2Read2Contig[sKey]:
				# FILE.write(sRead+"\t"+dPair2Read2Contig[sKey][sRead]+"\n")
		# else:
			# exit("ERROR-90: FATAL!")
	# FILE.close()
			
	# print("Pair mapped:",iPairOk)
	# print("Pair with ambigous mapping:",iPairAmbigous)
	# print("Pair with only one read mapped:",iPairHalfLost)
	# print("Pair without mapping",iPairLost)

	# return setLost,dContig2Read,dContig2Sample

def Parse(sPath):
	print("Parsing "+sPath)
	dContig2Read={}
	dContig2Sample={}
	dRead2Contig={}
	
	dPair2Reads={}

	iLineCounter=0

	for sNewLine in open(sPath):
		iLineCounter+=1
		if iLineCounter%50000==0:
			print("Reading line "+str(iLineCounter)+"...")
		if sNewLine[0]==REJECT_TAG:
			continue
		sLine=sNewLine.strip()
		tLine=sLine.split("\t")
		sMap=tLine[COL_MAPPING] #Match or not
		sId=tLine[COL_SEQID]	#Read id
		
		try:
			dContig2Read[sMap].add(sId)
		except KeyError:
			dContig2Read[sMap]=set([sId])
			
		try:
			dContig2Sample[sMap].add(sId.split("_")[-1])
		except KeyError:
			dContig2Sample[sMap]=set([sId.split("_")[-1]])
			
		try:
			dRead2Contig[sId].add(sMap)
		except KeyError:
			dRead2Contig[sId]=set([sMap])
			
		try:
			dPair2Reads[sId.split("_")[0]].add(sId)
		except KeyError:
			dPair2Reads[sId.split("_")[0]]=set([sId])
	
	setLost=set([])
	FILE=open(OUTPUT_AMBIGOUS_READ,"w")
	iReadOk=0
	iReadLost=0
	iReadAmbigous=0
	for sKey in dRead2Contig:
		if NO_MAPPING in dRead2Contig[sKey] and len(dRead2Contig[sKey])==1:
			iReadLost+=1
			setLost.add(sKey)
		elif NO_MAPPING in dRead2Contig[sKey] and len(dRead2Contig[sKey])!=1:
			exit("ERROR-198: {} have both NoMapping and Mapping".format(sKey))
		elif NO_MAPPING not in dRead2Contig[sKey] and len(dRead2Contig[sKey])==1:
			iReadOk+=1
		elif NO_MAPPING not in dRead2Contig[sKey] and len(dRead2Contig[sKey])!=1:
			iReadAmbigous+=1
			for sContig in sorted(dRead2Contig[sKey]):
				FILE.write(sKey+"\t"+sContig+"\n")
		else:
			exit("ERROR-90: FATAL!")
	FILE.close()
			
	print("Reads mapped:",iReadOk)
	print("Reads with ambigous mapping:",iReadAmbigous)
	print("Reads without mapping",iReadLost)
	
	iNoMapping=0
	iPartialMapping=0
	iPartialAmbigousMapping=0
	iSimilarMapping=0
	iDissimilarMapping=0
	iPartialSimilarAmbigousMapping=0
	iPartialDissimilarAmbigousMapping=0
	iSimilarAmbigousMapping=0
	iDissimilarAmbigousMapping=0
	iFuzzyAmbigousMapping=0
	
	for sPair in dPair2Reads:
		if len(dPair2Reads[sPair])==2:
			tPair=list(dPair2Reads[sPair])
			sPair1=tPair[0]
			sPair2=tPair[1]
			tTarget1=dRead2Contig[sPair1]
			tTarget2=dRead2Contig[sPair2]
			iCase1=0
			iCase2=0
			if NO_MAPPING in tTarget1:
				iCase1=1
			elif len(tTarget1)==1:
				iCase1=2
			elif len(tTarget1)>1:
				iCase1=3
			else:
				exit("ERROR-235: FATAL!")
			if NO_MAPPING in tTarget2:
				iCase2=1
			elif len(tTarget2)==1:
				iCase2=2
			elif len(tTarget2)>1:
				iCase2=3
			else:
				exit("ERROR-243: FATAL!")
			if iCase1==1 and iCase2==1:
				iNoMapping+=1
			elif 1 in [iCase1,iCase2] and 2 in [iCase1,iCase2]:
				iPartialMapping+=1
			elif 1 in [iCase1,iCase2] and 3 in [iCase1,iCase2]:
				iPartialAmbigousMapping+=1
			elif iCase1==2 and iCase2==2:
				if tTarget1[0]==tTarget2[0]:
					iSimilarMapping+=1
				else:
					iDissimilarMapping+=1
			elif 2 in [iCase1,iCase2] and 3 in [iCase1,iCase2]:
				if len(tTarget1)<len(tTarget2):
					sSingle=tTarget1[0]
					tMultiple=tTarget2
				else:
					sSingle=tTarget2[0]
					tMultiple=tTarget1
				if sSingle in tMultiple:
					iPartialSimilarAmbigousMapping+=1
				else:
					iPartialDissimilarAmbigousMapping+=1
			elif iCase1==3 and iCase2==3:
				tDelta=[X for X in tTarget1 if X not in tTarget2]
				if len(tDelta)==0:
					iSimilarAmbigousMapping+=1
				elif len(tDelta)==len(tTarget1):
					iDissimilarAmbigousMapping+=1
				else:
					iFuzzyAmbigousMapping+=1
			else:
				exit("ERROR-286: FATAL!")

	print("Pair fate:")
	print("NoMapping",iNoMapping,"(R1:None, R2:None)")
	print("PartialMapping",iPartialMapping,"(R1:None, R2:C1')")
	print("iPartialAmbigousMapping",iPartialAmbigousMapping,"(R1:None, R2:C1',C2')")
	print("iSimilarMapping",iSimilarMapping,"(R1:C1, R2:C1' ; C1==C1')")
	print("iDissimilarMapping",iDissimilarMapping,"(R1:C1, R2:C1' ; C1=/=C1')")
	print("iPartialSimilarAmbigousMapping",iPartialSimilarAmbigousMapping,"(R1:C1, R2:C1',C2' ; C1 in C1',C2')")
	print("iPartialDissimilarAmbigousMapping",iPartialSimilarAmbigousMapping,"(R1:C1, R2:C1',C2' ; C1 not in C1',C2')")
	print("iSimilarAmbigousMapping",iSimilarAmbigousMapping,"(R1:C1,C2, R2:C1',C2' ; C1,C2 == C1',C2')")
	print("iDissimilarAmbigousMapping",iDissimilarAmbigousMapping,"(R1:C1,C2, R2:C1',C2' ; C1,C2 =/= C1',C2')")
	print("iFuzzyAmbigousMapping",iFuzzyAmbigousMapping,"(R1:C1,C2, R2:C1',C2' ; C1 in C1',C2' and C2 not in C1',C2')")
	
	return setLost,dContig2Read,dContig2Sample

def WriteTab(dContig2Read,dContig2Sample):
	dContigId2NewId={}
	iId=0
	for sKey in dContig2Sample:
		if sKey==NO_MAPPING:
			continue
		iId+=1
		sName=CONTIG_BASENAME+"_"+str(iId)+"_("+str(len(dContig2Read[sKey]))+")"
		dContigId2NewId[sKey]=sName
	
	FILE=open(OUTPUT_FILENAME,"w")
	for sKey in dContig2Sample:
		if sKey==NO_MAPPING:
			continue
		tSample=dContig2Sample[sKey]
		tRead=dContig2Read[sKey]
		sName=dContigId2NewId[sKey]
		for sRead in sorted(tRead):
			FILE.write(sRead+"\t"+sName+"\t"+sKey+"\t"+"-".join(sorted(tSample))+"\n")
	FILE.close()
	
	return dContigId2NewId

def WriteContigFasta(dDict):
	FILE=open(OUTPUT_SPADES_CONTIG,"w")
	FILE2=open(OUTPUT_REJECTED_CONTIG,"w")
	bWriteIt=False
	for sNewLine in open(MEGAHIT_ASSEMBLY_OUTPUT):
		if ">" in sNewLine:
			sContigName=sNewLine[1:-1]
			sShortContigName=sContigName.split(" ")[0]
			try:
				sNewName=dDict[sShortContigName]
				bWriteIt=True
			except KeyError:
				bWriteIt=False
				print("Warning: Contig "+sContigName+" have no read reverse-mapped. Contig is removed.")
			if bWriteIt:
				FILE.write(">"+sNewName+"\n")
			else:
				FILE2.write(">"+sContigName+"\n")
		else:
			if bWriteIt:
				FILE.write(sNewLine)
			else:
				FILE2.write(sNewLine)
	FILE.close()
	FILE2.close()
		
def WriteLostData(dDict):
	setLost=set(dDict[NO_MAPPING])
	for sBase in dBASEtoTARGET:
		FILE=open(dBASEtoTARGET[sBase],"w")
		iLineCount=0
		bWriteIt=False
		for sNewLine in open(sBase):
			iLineCount+=1
			if iLineCount%4==1:
				bWriteIt=False
				sReadId=sNewLine[1:-1]
				if sReadId in setLost:
					bWriteIt=True
					FILE.write(sNewLine)
			elif bWriteIt:
					FILE.write(sNewLine)
		FILE.close()

########################################################################
#MAIN
if __name__ == "__main__":
	setLost,dContig2Read,dContig2Sample=Parse(sInput)
	print("len(setLost):",len(setLost))
	print("len(dContig2Read):",len(dContig2Read))
	print("len(dContig2Sample):",len(dContig2Sample))
	if len(dContig2Read)!=len(dContig2Sample):
		print("Spe dContig2Read:",[X for X in dContig2Read.keys() if X not in dContig2Sample])
		print("Spe dContig2Sample:",[X for X in dContig2Sample.keys() if X not in dContig2Read])
	dContig2Name=WriteTab(dContig2Read,dContig2Sample)
	print("len(dContig2Name):",len(dContig2Name))
	WriteContigFasta(dContig2Name)
	WriteLostData(dContig2Read)

########################################################################    
iTime2=time.time()
iDeltaTime=iTime2-iTime1
print("Script done: "+str(iDeltaTime))

