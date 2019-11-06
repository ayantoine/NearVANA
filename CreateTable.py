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
Create tsv bilan from Blast data

python CreateTable.py -t TASK -i JOBS -p PID
TASK: X or N
JOBS: jobs number
PID: Plaque Id
'''
########################################################################
#CONSTANT
HEADER_LIST=["Hit rank","Query Seq-Id","Sample","Read quantity","Sequence length","Location","Date","Host",
"Individual","Weight(mg)","Subject Seq-Id","Organism","SuperKingdom","Taxonomy","Hit definition","Identity",
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
parser.add_option("-t","--task", dest="task")
parser.add_option("-j","--jobs", dest="jobs")
parser.add_option("-p","--pid", dest="pid")
parser.add_option("-m","--meta", dest="meta")

(options, args) = parser.parse_args()

sTask=options.task
if not sTask:
	exit("Error : no task -t defined, process broken")

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

#Half-constant
BLAST_OUTPUT=sPID+"_Blast"+sTask+"_results.tab"
BLAST_FOLDER=sPID+"_Blast"+sTask
BLAST_INPUT=sPID+"_All.fa."+REPLACEME+".keeped"
BLAST_FILE=sPID+"_All.fa."+REPLACEME+".Blast"+sTask+"_2.tab"
TAXO_FILE=sPID+"_All.fa."+REPLACEME+".Blast"+sTask+"_2.tab.taxo"
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

def WriteData(FILE,dBlast,dTaxo,dContigs,dMetadata,dContent):
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
				
			for sSample in tSample:
				tLine=[sRank,sQuery,sSample,sReadQuantity,str(iQuerySize),
				dMetadata[sQuery]["Location"],dMetadata[sQuery]["Date"],
				dMetadata[sQuery]["Host"],dMetadata[sQuery]["Individuals"],
				dMetadata[sQuery]["Weight"],dBlast[sQuery][iRank]["SubjectId"],
				dTaxo[sQuery]["Organism"],dTaxo[sQuery]["Superkingdom"],
				dTaxo[sQuery]["Taxonomy"],dTaxo[sQuery]["Definition"],dBlast[sQuery][iRank]["Identity"],
				str(fCover),dBlast[sQuery][iRank]["Length"],dBlast[sQuery][iRank]["Mismatch"],
				dBlast[sQuery][iRank]["GapOpen"],dBlast[sQuery][iRank]["QueryStart"],
				dBlast[sQuery][iRank]["QueryEnd"],dBlast[sQuery][iRank]["SubjectStart"],
				dBlast[sQuery][iRank]["SubjectEnd"],dBlast[sQuery][iRank]["Evalue"],
				dBlast[sQuery][iRank]["BitScore"],dContent[sQuery]
				]
				FILE.write("\t".join(tLine)+"\n")
						
# HEADER_LIST=["Hit rank","Query Seq-Id","Sample","Read quantity","Sequence length","Location","Date","Host",
# "Individual","Weight(mg)","Subject Seq-Id","Organism","SuperKingdom","Taxonomy","Hit definition","Identity",
# "Query cover","Alignment length","Mismatches","Gap opening","Start alignment query","End alignment query",
# "Start alignment subject","End alignment subject","E-value","Bit score","Query sequences"]

########################################################################
#MAIN
if __name__ == "__main__":
	dMetadata=LoadMetadata(sMeta)
	dQuery2Content=LoadQuery(BLAST_INPUT)
	dContigs2Sample=LoadContigs(SHORTSPADES,dQuery2Size)
	dContigs2Sample=LoadContigs(SHORTFLASH,dQuery2Size,dContigs2Sample)
	FILE=open(BLAST_OUTPUT,"w")
	FILE.write(HEADER)
	for iIndex in range(1,iJobs+1):
		dTaxo=LoadTaxo(BLAST_FOLDER+"/"+TAXO_FILE.replace(REPLACEME,str(iIndex)))
		dBlast=LoadBlast(BLAST_FOLDER+"/"+BLAST_FILE.replace(REPLACEME,str(iIndex)))
		WriteData(FILE,dBlast,dTaxo,dContigs2Sample,dMetadata,dQuery2Content)
	FILE.close()
	
########################################################################    
iTime2=time.time()
iDeltaTime=iTime2-iTime1
print("Script done: "+str(iDeltaTime))

