# coding: utf-8
"""Python3.6"""
# compatibility: python2.7, python2.6

import time
from optparse import OptionParser

sCurrentVersionScript="v1"
iTime1=time.time()
########################################################################
'''
V1-2019/11/19
Count and store numbrer of read by sample from different files

python Compute.py -p PID
PID: Plaque Id
'''
########################################################################
#CONSTANT
DEMULTIPLEXING_SUFFIX="_Demultiplexing_Hyper_Distribution.tsv"
ALLSEQ_SUFFIX="_All.fa"
BLASTXKEEPED_SUFFIX="_BlastX.keeped.fa"
BLASTNKEEPED_SUFFIX="_BlastN.keeped.fa"
FLASHREVERSE_SUFFIX="_All.FLASH_reverseAssembly.tsv"
SPADESREVERSE_SUFFIX="_All.SPAdes_reverseAssembly.tsv"

STATISTIC_SUFFIX="_Statistic.tsv"

OUTPUT_SUFFIX="_Statistics.tsv"

CONTIG="Contigs"

########################################################################
#Options
if __name__ == "__main__":
	parser = OptionParser()
	parser.add_option("-p","--pid", dest="pid")

	(options, args) = parser.parse_args()

	sPid=options.pid
	if not sPid:
		exit("Error : no pid -p defined, process broken")

########################################################################
#Function 	
def LoadDecomplexing(sPid,sPath):
	dResult={"Demultiplexing":{},"Sample":{}}
	bHeader=True
	print("Loading "+sPath)
	for sNewLine in open(sPath):
		if bHeader:
			bHeader=False
			continue
		sLine=sNewLine.strip()
		tLine=sLine.split("\t")
		sSample=tLine[0].replace(sPid+"_","")
		sQuantity=tLine[1]
		dResult["Sample"][sSample]=1
		dResult["Demultiplexing"][sSample]=sQuantity
	return dResult
		
def LoadReverseMapping(sPath,dDict={}):
	print("Loading "+sPath)
	for sNewLine in open(sPath):
		sLine=sNewLine.strip()
		tLine=sLine.split("\t")
		sRead=tLine[0]
		sSample=sRead.split("_")[-1]
		sContig=tLine[1]
		try:
			oCrash=dDict[sContig]
		except KeyError:
			dDict[sContig]={}
		try:
			dDict[sContig][sSample]+=1
		except KeyError:
			dDict[sContig][sSample]=1
	return dDict

def LoadAllSeq(dResult,dDict,sPath):
	print("Loading "+sPath)
	dResult["AllSeq"]={}
	for sNewLine in open(sPath):
		if ">"==sNewLine[0]:
			sLine=sNewLine.strip()
			if CONTIG in sLine:
				sContig=sLine[1:]
				try:
					oCrash=dDict[sContig]
				except KeyError:
					print("{} not found in data".format(sContig))
				for sSample in dDict[sContig]:
					try:
						dResult["AllSeq"][sSample]+=dDict[sContig][sSample]
					except KeyError:
						dResult["AllSeq"][sSample]=dDict[sContig][sSample]
					try:
						oCrash=dResult["Sample"][sSample]
					except KeyError:
						dResult["Sample"][sSample]=1
			else:
				sSample=sLine.split("_")[-1]
				try:
					dResult["AllSeq"][sSample]+=1
				except KeyError:
					dResult["AllSeq"][sSample]=1
				try:
					oCrash=dResult["Sample"][sSample]
				except KeyError:
					dResult["Sample"][sSample]=1
	return dResult
	
def LoadBlastNKeeped(dResult,dDict,sPath):
	print("Loading "+sPath)
	dResult["BlastNKeeped"]={}
	print("oco")
	for sNewLine in open(sPath):
	# for sNewLine in open(sPath):
		print(sNewLine)
		if ">"==sNewLine[0]:
			print(sNewLine)
			sLine=sNewLine.strip()
			if CONTIG in sLine:
				sContig=sLine[1:]
				print(sContig)
				try:
					oCrash=dDict[sContig]
				except KeyError:
					print("{} not found in data".format(sContig))
				for sSample in dDict[sContig]:
					try:
						oCrash=dResult["Sample"][sSample]
					except KeyError:
						dResult["Sample"][sSample]=1
					try:
						dResult["BlastNKeeped"][sSample]+=dDict[sContig][sSample]
					except KeyError:
						dResult["BlastNKeeped"][sSample]=dDict[sContig][sSample]
			else:
				print("Not a contig")
				sSample=sLine.split("_")[-1]
				try:
					oCrash=dResult["Sample"][sSample]
				except KeyError:
					dResult["Sample"][sSample]=1
				try:
					dResult["BlastNKeeped"][sSample]+=1
				except KeyError:
					dResult["BlastNKeeped"][sSample]=1
	return dResult
	
def LoadBlastXKeeped(dResult,dDict,sPath):
	print("Loading "+sPath)
	dResult["BlastXKeeped"]={}
	for sNewLine in open(sPath):
		print(sNewLine)
		if ">"==sNewLine[0]:
			sLine=sNewLine.strip()
			if CONTIG in sLine:
				sContig=sLine[1:]
				try:
					oCrash=dDict[sContig]
				except KeyError:
					print("{} not found in data".format(sContig))
				for sSample in dDict[sContig]:
					try:
						oCrash=dResult["Sample"][sSample]
					except KeyError:
						dResult["Sample"][sSample]=1
					try:
						dResult["BlastXKeeped"][sSample]+=dDict[sContig][sSample]
					except KeyError:
						dResult["BlastXKeeped"][sSample]=dDict[sContig][sSample]
			else:
				sSample=sLine.split("_")[-1]
				try:
					oCrash=dResult["Sample"][sSample]
				except KeyError:
					dResult["Sample"][sSample]=1
				try:
					dResult["BlastXKeeped"][sSample]+=1
				except KeyError:
					dResult["BlastXKeeped"][sSample]=1
	return dResult
		
def WriteStat(dResult,sPath):
	FILE=open(sPath,"w")
	tCategorie=["Sample","Demultiplexing","AllSeq","BlastNKeeped","BlastXKeeped"]
	FILE.write("\t".join(tCategorie)+"\n")
	for sSample in sorted(dResult[tCategorie[0]]):
		FILE.write(sSample)
		for sCategorie in tCategorie[1:]:
			FILE.write("\t")
			try:
				iValue=dResult[sCategorie][sSample]
			except KeyError:
				iValue=0
			FILE.write(str(iValue))
		FILE.write("\n")
	FILE.close()

########################################################################
#MAIN
if __name__ == "__main__":
	dReverseMapping=LoadReverseMapping(sPid+SPADESREVERSE_SUFFIX)
	dReverseMapping=LoadReverseMapping(sPid+FLASHREVERSE_SUFFIX,dReverseMapping)
	dData=LoadDecomplexing(sPid,sPid+DEMULTIPLEXING_SUFFIX)
	dData=LoadAllSeq(dData,dReverseMapping,sPid+ALLSEQ_SUFFIX)
	dData=LoadBlastNKeeped(dData,dReverseMapping,sPid+BLASTNKEEPED_SUFFIX)
	dData=LoadBlastXKeeped(dData,dReverseMapping,sPid+BLASTXKEEPED_SUFFIX)
	WriteStat(dData,sPid+OUTPUT_SUFFIX)
	
########################################################################    
iTime2=time.time()
iDeltaTime=iTime2-iTime1
print("Script done: "+str(iDeltaTime))

