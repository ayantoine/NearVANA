# coding: utf-8
"""Python3.6"""
# compatibility: python2.7, python2.6

import time
from optparse import OptionParser

sCurrentVersionScript="v1"
iTime1=time.time()
########################################################################
'''
V1-2019/11/06
Create tsv bilan from BlastX and BlastN data

python CreateTable.py -i JOBS -p PID
JOBS: jobs number
PID: Plaque Id
'''
########################################################################
#CONSTANT
HEADER_LIST=["Hit rank","Query Seq-Id","Sample","Read quantity","Sequence length","Location","Date","Host",
"Individual","Weight(mg)","Subject Seq-Id","Organism","SuperKingdom","Taxonomy","Hit definition","Identification","% Fragment","Identity",
"Query cover","Alignment length","Mismatches","Gap opening","Start alignment query","End alignment query",
"Start alignment subject","End alignment subject","E-value","Bit score","Query sequences"]

HEADER="\t".join(HEADER_LIST)+"\n"

BEST_HIT="Best hit"
HIT="."

REPLACEME="REPLACE-ME"
SAMPLE_SEPARATOR="-"
CONTIG="Contig"

META_SAMPLECOL=0
META_HOSTCOL=1
META_LOCATIONCOL=3
META_DATECOL=6
META_INDIVIDUALSCOL=7
META_WEIGHTCOL=8

BLAST_QUERYIDCOl=0
BLAST_SUBJECTIDCOl=1
BLAST_IDENTITYCOl=2
BLAST_LENGTHCOl=3
BLAST_MISMATCHCOl=4
BLAST_GAPOPENCOl=5
BLAST_QUERYSTARTCOl=6
BLAST_QUERYENDCOl=7
BLAST_SUBJECTSTARTCOl=8
BLAST_SUBJECTENDCOl=9
BLAST_EVALUECOl=10
BLAST_BITSCORECOl=11

TAXO_ACCCOL=0
TAXO_ORGANISMCOL=1
TAXO_SUPKINGDOMCOL=2
TAXO_LINEAGECOL=3
TAXO_DEFCOL=4

########################################################################
#Options
parser = OptionParser()
parser.add_option("-j","--jobs", dest="jobs")
parser.add_option("-p","--pid", dest="pid")
parser.add_option("-m","--meta", dest="meta")
parser.add_option("-l","--length", dest="length")

(options, args) = parser.parse_args()

sJobs=options.jobs
if not sJobs:
	exit("Error : no jobs -j defined, process broken")
try:
	iJobs=int(sJobs)
except KeyError:
	exit("Error : jobs -j must be an integer, process broken")
	
sPID=options.pid
if not sPID:
	sys.exit("Error : no pid -p defined, process broken")

sMeta=options.meta
if not sMeta:
	sys.exit("Error : no meta -m defined, process broken")


sLengthFile=options.length
if not sLengthFile:
	exit("Error : no length -l defined, process broken")

#Half-constant
BLAST_INPUT=sPID+"_All.fa."+REPLACEME+".keeped"

BLASTALL_OUTPUT=sPID+"_BlastAll_results.tab"

BLASTN_FOLDER=sPID+"_BlastN"
BLASTN_FILE=sPID+"_All.fa."+REPLACEME+".BlastN_2.tab"
TAXON_FILE=sPID+"_All.fa."+REPLACEME+".BlastN_2.tab.taxo"

BLASTX_FOLDER=sPID+"_BlastX"
BLASTX_FILE=sPID+"_All.fa."+REPLACEME+".BlastX_2.tab"
TAXOX_FILE=sPID+"_All.fa."+REPLACEME+".BlastX_2.tab.taxo"

SHORTSPADES=sPID+"_All.SPAdes.contigs2sample.tsv"
SHORTFLASH=sPID+"_All.FLASH.contigs2sample.tsv"

########################################################################
#Function 	
def LoadMetadata(sFile):
	dDict={}
	for sNewLine in open(sFile):
		sLine=sNewLine.strip()
		tLine=sLine.split()
		dDict[tLine[META_SAMPLECOL]]={
			"Host":tLine[META_HOSTCOL],
			"Location":tLine[META_LOCATIONCOL],
			"Date":tLine[META_DATECOL],
			"Individuals":tLine[META_INDIVIDUALSCOL],
			"Weight":tLine[META_WEIGHTCOL]
			}
	return dDict

def LoadContigs(sFile,dRef,dDict={}):
	for sNewLine in open(sFile):
		sLine=sNewLine.strip()
		tLine=sLine.split()
		sName=tLine[0]
		try:
			oCrash=dRef[sName]
			dDict[sName]=tLine[1].split(SAMPLE_SEPARATOR)
		except KeyError:
			continue
	return dDict

def LoadTaxo(sFile):
	dDict={}
	for sNewLine in open(sFile):
		sLine=sNewLine.strip()
		tLine=sLine.split()
		dDict[tLine[TAXO_ACCCOL]]={
			"Organism":tLine[TAXO_ORGANISMCOL],
			"Superkingdom":tLine[TAXO_SUPKINGDOMCOL],
			"Lineage":tLine[TAXO_LINEAGECOL],
			"Definition":tLine[TAXO_DEFCOL]
			}
	return dDict

def LoadBlast(sFile):
	dDict={}
	for sNewLine in open(sFile):
		sLine=sNewLine.strip()
		tLine=sLine.split()
		sQueryId=tLine[BLAST_QUERYIDCOl]
		try:
			oCrash=dDict[sQueryId]
		except KeyError:
			dDict[sQueryId]={}
		dDict[sQueryId][len(dDict[sQueryId])+1]={
			"SubjectId":tLine[BLAST_SUBJECTIDCOl],
			"Identity":tLine[BLAST_IDENTITYCOl],
			"Length":tLine[BLAST_LENGTHCOl],
			"Mismatch":tLine[BLAST_MISMATCHCOl],
			"GapOpen":tLine[BLAST_GAPOPENCOl],
			"QueryStart":tLine[BLAST_QUERYSTARTCOl],
			"QueryEnd":tLine[BLAST_QUERYENDCOl],
			"SubjectStart":tLine[BLAST_SUBJECTSTARTCOl],
			"SubjectEnd":tLine[BLAST_SUBJECTENDCOl],
			"Evalue":tLine[BLAST_EVALUECOl],
			"BitScore":tLine[BLAST_BITSCORECOl]
			}
	return dDict

def LoadQuery(sFile):
	dDict={}
	sSeqContent=""
	sSeqName=""
	for sNewLine in open(sFile):
		if ">"==sNewLine[0]:
			if sSeqName!="":
				dDict[sSeqName]=sSeqContent
			sSeqName=sNewLine[1:-1]
			sSeqContent=""
		else:
			sSeqContent+=sNewLine[:-1] #:-1, for \n char
	if ">"==sNewLine[0]:
		if sSeqName!="":
			dDict[sSeqName]=sSeqContent
	return dDict

def WriteData(FILE,dBlast,dTaxo,dContigs,dMetadata,dContent,dLength):
	for sQuery in dBlast:
		for iRank in dBlast[sQuery]:
			
			if iRank==1:
				sRank=BEST_HIT
			else:
				sRank=HIT
			
			iQuerySize=len(dContent[sQuery])
			iCoverSize=int(dBlast[sQuery][iRank][QueryEnd])-int(dBlast[sQuery][iRank][QueryStart])
			fCover=round(float(iCoverSize)/iQuerySize*100,2)
			
			try:
				tSample=dContigs[sQuery]
			except KeyError:
				tSample=[sQuery.split("_")[-1]]
			
			sReadQuantity=sQuery.split("(")[-1].split(")")[0]
			sSubjectId=dBlast[sQuery][iRank]["SubjectId"]
			
			sTaxo=dTaxo[sSubjectId]["Lineage"]
			tTaxo=sTaxo.replace("; ",";").split(";")
			iMinSize=0
			for oTaxo in tTaxo:
				try:
					iMinSize=dLength[oTaxo]
				except KeyError:
					continue
			if iMinSize==0:
				fFragment="N/A"
			else:
				fFragment=round(float(iQuerySize)/iMinSize*100,2)
				
			for sSample in tSample:
				tLine=[sRank,sQuery,sSample,sReadQuantity,str(iQuerySize),
				dMetadata[sSample]["Location"],dMetadata[sSample]["Date"],
				dMetadata[sSample]["Host"],dMetadata[sSample]["Individuals"],
				dMetadata[sSample]["Weight"],sSubjectId,
				dTaxo[sSubjectId]["Organism"],dTaxo[sSubjectId]["Superkingdom"],
				dTaxo[sSubjectId]["Lineage"],dTaxo[sSubjectId]["Definition"],dBlast[sQuery][iRank]["Confidence"],
				str(fFragment),dBlast[sQuery][iRank]["Identity"],
				str(fCover),dBlast[sQuery][iRank]["Length"],dBlast[sQuery][iRank]["Mismatch"],
				dBlast[sQuery][iRank]["GapOpen"],dBlast[sQuery][iRank]["QueryStart"],
				dBlast[sQuery][iRank]["QueryEnd"],dBlast[sQuery][iRank]["SubjectStart"],
				dBlast[sQuery][iRank]["SubjectEnd"],dBlast[sQuery][iRank]["Evalue"],
				dBlast[sQuery][iRank]["BitScore"],dContent[sQuery]
				]
				FILE.write("\t".join(tLine)+"\n")
						
# HEADER_LIST=["Hit rank","Query Seq-Id","Sample","Read quantity","Sequence length","Location","Date","Host",
# "Individual","Weight(mg)","Subject Seq-Id","Organism","SuperKingdom","Taxonomy","Hit definition","Identification","% Fragment","Identity",
# "Query cover","Alignment length","Mismatches","Gap opening","Start alignment query","End alignment query",
# "Start alignment subject","End alignment subject","E-value","Bit score","Query sequences"]

def FusionTaxoDict(dDict1,dDict2):
	dDict={}
	for sKey in dDict1:
		dDict[sKey]=dDict1[sKey]
	for sKey in dDict2:
		dDict[sKey]=dDict2[sKey]
	return dDict
	
def FusionBlastDict(dBlastN,dBlastX):
	dDict={}
	tKey=list(set(list(dBlastN.keys())+list(dBlastX.keys())))
	
	for sKey in tKey:
		bX=True
		bN=True
		try:
			dX=dBlastX[sKey]
		except KeyError:
			bX=False
		try:
			dN=dBlastN[sKey]
		except KeyError:
			bN=False
			
		if bX and not bN:
			#Identified only with BlastX
			for iRank in dBlastX[sKey]:
				dDict[sKey][iRank]=dBlastX[sKey][iRank]
				dDict[sKey][iRank]["Confidence"]="Single (X)"
	
		elif not bX and bN:
			#Identified only with BlastN
			for iRank in dBlastN[sKey]:
				dDict[sKey][iRank]=dBlastN[sKey][iRank]
				dDict[sKey][iRank]["Confidence"]="Single (N)"
				
		else:
			dDict[sKey]={}
			#Identified with both BlastX and BlastN
			dSubjectId2MaxBitScore={}
			dSubjectId2Confidence={}
			dSubjectId2Coord={}
			for sKey in dBlastN:
				for iRank in dBlastN[sKey]:
					dSubjectId2MaxBitScore[dBlastN[sKey][iRank]["SubjectId"]]=float(dBlastN[sKey][iRank]["BitScore"])
					dSubjectId2Confidence[dBlastN[sKey][iRank]["SubjectId"]]="Single (N)"
					dSubjectId2Coord[dBlastN[sKey][iRank]["SubjectId"]]=("N",iRank)
			for sKey in dBlastX:
				for iRank in dBlastX[sKey]:
					sSubjectId=dBlastX[sKey][iRank]["SubjectId"]
					fBitScore=float(dBlastX[sKey][iRank]["BitScore"])
					if sSubjectId in dSubjectId2MaxBitScore:
						dSubjectId2Confidence[sSubjectId]="Double"
						# dSubjectId2MaxBitScore[sSubjectId]=max(dSubjectId2MaxBitScore[sSubjectId],fBitScore)
						if fBitScore>dSubjectId2MaxBitScore[sSubjectId]:
							dSubjectId2MaxBitScore[sSubjectId]=fBitScore
							dSubjectId2Coord[dBlastN[sKey][iRank]["SubjectId"]]=("X",iRank)
					else:
						dSubjectId2MaxBitScore[dBlastX[sKey][iRank]["SubjectId"]]=float(dBlastX[sKey][iRank]["BitScore"])
						dSubjectId2Confidence[dBlastX[sKey][iRank]["SubjectId"]]="Single (X)"
						dSubjectId2Coord[dBlastN[sKey][iRank]["SubjectId"]]=("X",iRank)
			dBitScore2SubjectId={}
			for sKey in dSubjectId2MaxBitScore:
				dBitScore2SubjectId[dSubjectId2MaxBitScore[sKey]]=sKey
			for fBitScore in sorted(dBitScore2SubjectId,reverse=True):
				sSubjectId=dBitScore2SubjectId[fBitScore]
				sConfidence=dSubjectId2Confidence[sSubjectId]
				dbCoord=dSubjectId2Coord[sSubjectId]
				iRank=dbCoord[1]
				if dbCoord[0]=="X":
					iNewRank=len(dDict[sKey])+1
					dDict[sKey][iNewRank]=dBlastX[sSubjectId][iRank]
					dDict[sKey][iNewRank]["Confidence"]=dSubjectId2Confidence[sSubjectId]
				else:
					iNewRank=len(dDict[sKey])+1
					dDict[sKey][iNewRank]=dBlastN[sSubjectId][iRank]
					dDict[sKey][iNewRank]["Confidence"]=dSubjectId2Confidence[sSubjectId]
					
	return dDict

def LoadLength(sFile):
	print("Loading file "+str(sFile))
	dDict={}
	for sNewLine in open(sFile):
		sLine=sNewLine.strip()
		tLine=sLine.split()
		sFamily=tLine[0]
		iMinSize=int(tLine[1])
		dDict[sFamily]=iMinSize
	return dDict					

########################################################################
#MAIN
if __name__ == "__main__":
	dMetadata=LoadMetadata(sMeta)
	dLength=LoadLength(sLengthFile)
	FILE=open(BLAST_OUTPUT,"w")
	FILE.write(HEADER)
	for iIndex in range(1,iJobs+1):
		dQuery2Content=LoadQuery(BLAST_FOLDER+"/"+BLAST_INPUT.replace(REPLACEME,str(iIndex)))
		dContigs2Sample=LoadContigs(SHORTSPADES,dQuery2Content)
		dContigs2Sample=LoadContigs(SHORTFLASH,dQuery2Content,dContigs2Sample)
		dTaxoN=LoadTaxo(BLASTN_FOLDER+"/"+TAXON_FILE.replace(REPLACEME,str(iIndex)))
		dTaxoX=LoadTaxo(BLASTX_FOLDER+"/"+TAXOX_FILE.replace(REPLACEME,str(iIndex)))
		dBlastN=LoadBlast(BLASTN_FOLDER+"/"+BLASTN_FILE.replace(REPLACEME,str(iIndex)))
		dBlastX=LoadBlast(BLASTX_FOLDER+"/"+BLASTX_FILE.replace(REPLACEME,str(iIndex)))
		dTaxo=FusionTaxoDict(dTaxoN,dTaxoX)
		dBlast=FusionBlastDict(dBlastN,dBlastX)
		WriteData(FILE,dBlast,dTaxo,dContigs2Sample,dMetadata,dQuery2Content,dLength)
	FILE.close()
	
########################################################################    
iTime2=time.time()
iDeltaTime=iTime2-iTime1
print("Script done: "+str(iDeltaTime))

