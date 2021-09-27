# coding: utf-8
"""Python3.6"""
# compatibility: python2.7, python2.6

import time
from optparse import OptionParser

sCurrentVersionScript="v1"
iTime1=time.time()
########################################################################
'''
V4-2021/09/24


V3-2021/04/12
Fix: No more able to catch metadata file
V2-2020/02/14
Adapt to Diamond and Demultiplexing
V1-2019/11/06
Create tsv bilan from Blast data

python CreateTable.py -d DATA -i JOBS -p PID
DATA: Data argfile
JOBS: jobs number
PID: Plaque Id
'''
########################################################################
#CONSTANT
HEADER_LIST=["Hit rank","Query Seq-Id","Sample","Read quantity","Sequence length","Location","Date","Host",
"Individual","Weight(mg)","Subject Seq-Id","Organism","SuperKingdom","Taxonomy","Hit definition","% Fragment","Identity",
"Query cover","Alignment length","Mismatches","Gap opening","Start alignment query","End alignment query",
"Start alignment subject","End alignment subject","E-value","Bit score","Query sequences"]
HEADER="\t".join(HEADER_LIST)+"\n"

BEST_HIT="Best hit"
HIT="."

REPLACEME="REPLACE-ME"
SAMPLE_SEPARATOR="-"
CONTIG="Contig"
TABULATION="\t"

UNASSIGNED_READS="UnassignedReads"

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

DEFAULT="."

EMPTY=""
DIESE="#"
EQUAL="="
OPEN_PARENTHESIS="("
CLOSE_PARENTHESIS=")"
SPACE=" "
PIPE="|"

########################################################################
#Options
if __name__ == "__main__":
	parser = OptionParser()
	parser.add_option("-j","--jobs", dest="jobs")
	parser.add_option("-p","--pid", dest="pid")
	parser.add_option("-l","--length", dest="length")
	parser.add_option("-d","--datafile", dest="datafile")
	parser.add_option("-t","--task", dest="task")

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
		exit("Error : no pid -p defined, process broken")

	sLengthFile=options.length
	if not sLengthFile:
		exit("Error : no length -l defined, process broken")
	
	sDataFile=options.datafile
	if not sDataFile:
		exit("Error : no datafile -d defined, process broken")
		
	sTask=options.task
	if not sTask:
		exit("Error : no task -t defined, process broken")

	if sTask=="D":
		#Half-constant
		BLAST_OUTPUT=sPID+"_Blast"+sTask+"_results.tab"
		BLAST_FOLDER=sPID+"_Blast"+sTask
		BLAST_INPUT=sPID+"_All.fa."+REPLACEME+".keeped"
		BLAST_FILE=sPID+"_All.fa."+REPLACEME+".Diamond_2.tab"
		TAXO_FILE=sPID+"_All.fa."+REPLACEME+".Diamond_2.tab.taxo"
	else:
		#Half-constant
		BLAST_OUTPUT=sPID+"_Blast"+sTask+"_results.tab"
		BLAST_FOLDER=sPID+"_Blast"+sTask
		BLAST_INPUT=sPID+"_All.fa."+REPLACEME+".keeped"
		BLAST_FILE=sPID+"_All.fa."+REPLACEME+".Blast"+sTask+"_2.tab"
		TAXO_FILE=sPID+"_All.fa."+REPLACEME+".Blast"+sTask+"_2.tab.taxo"
	REVERSE_ASSEMBLY=sPID+"_All.Megahit_reverseAssembly.tsv"

########################################################################
#Function 	
def LoadData(sFile):
	print("Loading file "+str(sFile))
	dDict={}
	bFirst=True
	for sNewLine in open(sFile):
		if sNewLine[0]==DIESE:
			continue
		if bFirst:
			bFirst=False
			continue
		sLine=sNewLine.strip()
		if sLine==EMPTY:
			continue
		tLine=sLine.split(EQUAL)
		sPlateId=tLine[0]
		sListOfData=tLine[1].replace(OPEN_PARENTHESIS,EMPTY).replace(CLOSE_PARENTHESIS,EMPTY)
		tListOfData=sListOfData.split(SPACE)
		sMeta=tListOfData[-1]
		dDict[sPlateId]=sMeta
	return dDict

def LoadMetadata(dData):
	dDict={}
	for sPlateId in dData:
		dDict[sPlateId]={}
		sFile=dData[sPlateId]
		print("Loading file "+str(sFile))
		for sNewLine in open(sFile):
			sLine=sNewLine.replace("\n","")
			if len(sLine)==0:
				continue
			tLine=sLine.split("\t")
						
			sHost=DEFAULT
			if tLine[META_HOSTCOL]!="":
				sHost=tLine[META_HOSTCOL]
			sLocation=DEFAULT
			if tLine[META_LOCATIONCOL]!="":
				sLocation=tLine[META_LOCATIONCOL]
			sDate=DEFAULT
			if tLine[META_DATECOL]!="":
				sDate=tLine[META_DATECOL]
			sIndividuals=DEFAULT
			if tLine[META_INDIVIDUALSCOL]!="":
				sIndividuals=tLine[META_INDIVIDUALSCOL]
			sWeight=DEFAULT
			if tLine[META_WEIGHTCOL]!="":
				sWeight=tLine[META_WEIGHTCOL]
			
			dDict[sPlateId][tLine[META_SAMPLECOL]]={
				"Host":sHost,
				"Location":sLocation,
				"Date":sDate,
				"Individuals":sIndividuals,
				"Weight":sWeight
				}
				
	return dDict

def LoadContigsAndQuantity(sFileReverse,dRef,dDict={}):
	print("Loading targeting reverse query from "+str(sFileReverse))
	for sNewLine in open(sFileReverse):
		sLine=sNewLine.strip()
		tLine=sLine.split(TABULATION)
		sReadId=tLine[0]
		sContigId=tLine[1]
		try:
			oCrash=dRef[sContigId]
			sSampleId=sReadId.split(TABULATION)[-1]
			if sContigId not in dDict:
				dDict[sContigId]={}
			try:
				dDict[sContigId][sSampleId]+=1
			except KeyError:
				dDict[sContigId][sSampleId]=1
		except KeyError:
			continue
	return dDict

def LoadTaxo(sFile):
	print("Loading file "+str(sFile))
	dDict={}
	for sNewLine in open(sFile):
		sLine=sNewLine.strip()
		tLine=sLine.split("\t")
		if len(tLine)!=5:
			continue
		
		# print(sFile)
		# print(tLine)
		sOrganism=DEFAULT
		if tLine[TAXO_ORGANISMCOL]!="":
			sOrganism=tLine[TAXO_ORGANISMCOL]
		sSuperkingdom=DEFAULT
		if tLine[TAXO_SUPKINGDOMCOL]!="":
			sSuperkingdom=tLine[TAXO_SUPKINGDOMCOL]
		sLineage=DEFAULT
		if tLine[TAXO_LINEAGECOL]!="":
			sLineage=tLine[TAXO_LINEAGECOL]
		sDefinition=DEFAULT
		if tLine[TAXO_DEFCOL]!="":
			sDefinition=tLine[TAXO_DEFCOL]
		
		
		dDict[tLine[TAXO_ACCCOL]]={
			"Organism":sOrganism,
			"Superkingdom":sSuperkingdom,
			"Lineage":sLineage,
			"Definition":sDefinition
			}
			
	return dDict

def LoadBlast(sFile):
	print("Loading file "+str(sFile))
	dDict={}
	try:
		for sNewLine in open(sFile):
			sLine=sNewLine.strip()
			tLine=sLine.split("\t")
			sQueryId=tLine[BLAST_QUERYIDCOl]
			try:
				oCrash=dDict[sQueryId]
			except KeyError:
				dDict[sQueryId]={}
			
			sSubjectId=DEFAULT
			sSubjectId=tLine[BLAST_SUBJECTIDCOl]
			if PIPE in sSubjectId:
				if PIPE==sSubjectId[-1]:
					sSubjectId=sSubjectId[:-1]
				sSubjectId=sSubjectId.split(PIPE)[-1]
			if sSubjectId==EMPTY:
				continue
			sIdentity=DEFAULT
			if tLine[BLAST_IDENTITYCOl]!="":
				sIdentity=tLine[BLAST_IDENTITYCOl]
			sLength=DEFAULT
			if tLine[BLAST_LENGTHCOl]!="":
				sLength=tLine[BLAST_LENGTHCOl]
			sMismatch=DEFAULT
			if tLine[BLAST_MISMATCHCOl]!="":
				sMismatch=tLine[BLAST_MISMATCHCOl]
			sGapOpen=DEFAULT
			if tLine[BLAST_GAPOPENCOl]!="":
				sGapOpen=tLine[BLAST_GAPOPENCOl]
			sQueryStart=DEFAULT
			if tLine[BLAST_QUERYSTARTCOl]!="":
				sQueryStart=tLine[BLAST_QUERYSTARTCOl]
			sQueryEnd=DEFAULT
			if tLine[BLAST_QUERYENDCOl]!="":
				sQueryEnd=tLine[BLAST_QUERYENDCOl]
			sSubjectStart=DEFAULT
			if tLine[BLAST_SUBJECTSTARTCOl]!="":
				sSubjectStart=tLine[BLAST_SUBJECTSTARTCOl]
			sSubjectEnd=DEFAULT
			if tLine[BLAST_SUBJECTENDCOl]!="":
				sSubjectEnd=tLine[BLAST_SUBJECTENDCOl]
			sEvalue=DEFAULT
			if tLine[BLAST_EVALUECOl]!="":
				sEvalue=tLine[BLAST_EVALUECOl]
			sBitScore=DEFAULT
			if tLine[BLAST_BITSCORECOl]!="":
				sBitScore=tLine[BLAST_BITSCORECOl]
						
			dDict[sQueryId][len(dDict[sQueryId])+1]={
				"SubjectId":sSubjectId,
				"Identity":sIdentity,
				"Length":sLength,
				"Mismatch":sMismatch,
				"GapOpen":sGapOpen,
				"QueryStart":sQueryStart,
				"QueryEnd":sQueryEnd,
				"SubjectStart":sSubjectStart,
				"SubjectEnd":sSubjectEnd,
				"Evalue":sEvalue,
				"BitScore":sBitScore
				}
	except FileNotFoundError:
		pass
	return dDict

def LoadQuery(sFile):
	print("Loading file "+str(sFile))
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
	if sSeqName!="":
		dDict[sSeqName]=sSeqContent
	return dDict

def WriteData(FILE,dBlast,dTaxo,dContigs,dMetadata,dContent,dLength):
	for sQuery in dBlast:
		# print("sQuery",sQuery)
		for iRank in dBlast[sQuery]:
			# print("iRank",iRank)
			if iRank==1:
				sRank=BEST_HIT
			else:
				sRank=HIT
			
			iQuerySize=len(dContent[sQuery])
			iCoverSize=int(dBlast[sQuery][iRank]["QueryEnd"])-int(dBlast[sQuery][iRank]["QueryStart"])
			fCover=round(float(iCoverSize)/iQuerySize*100,2)
			
			tGlobalSample=sorted(list(dContigs[sQuery].keys()))

						
			# if CONTIG in sQuery:
				# sReadQuantity=sQuery.split("(")[-1].split(")")[0]
			# else:
				# sReadQuantity="1"
			sSubjectId=dBlast[sQuery][iRank]["SubjectId"]
			
			# print("sReadQuantity",sReadQuantity)
			# print("sSubjectId",sSubjectId)
			
			try:
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
					
				sOrganism=dTaxo[sSubjectId]["Organism"]
				sSuperkingdom=dTaxo[sSubjectId]["Superkingdom"]
				sLineage=dTaxo[sSubjectId]["Lineage"]
				sDefinition=dTaxo[sSubjectId]["Definition"]

			except KeyError:
				sTaxo="unknown"
				fFragment="unknown"
				sOrganism="unknown"
				sSuperkingdom="unknown"
				sLineage="unknown"
				sDefinition="unknown"
				
			for sGlobalSample in tGlobalSample:
				# print("sGlobalSample",sGlobalSample)
				if sGlobalSample!=UNASSIGNED_READS:
					for sPlateId in dMetadata:
						# print("sPlateId",sPlateId)
						if sPlateId in sGlobalSample:
							break
					
					sReadQuantity=dContigs[sGlobalSample]
					sSampleId=sGlobalSample.replace(sPlateId,EMPTY)
					tLine=[sRank,sQuery,sGlobalSample,sReadQuantity,str(iQuerySize),
					dMetadata[sPlateId][sSampleId]["Location"],dMetadata[sPlateId][sSampleId]["Date"],
					dMetadata[sPlateId][sSampleId]["Host"],dMetadata[sPlateId][sSampleId]["Individuals"],
					dMetadata[sPlateId][sSampleId]["Weight"],sSubjectId,
					sOrganism,sSuperkingdom,
					sLineage,sDefinition,str(fFragment),
					dBlast[sQuery][iRank]["Identity"],
					str(fCover),dBlast[sQuery][iRank]["Length"],dBlast[sQuery][iRank]["Mismatch"],
					dBlast[sQuery][iRank]["GapOpen"],dBlast[sQuery][iRank]["QueryStart"],
					dBlast[sQuery][iRank]["QueryEnd"],dBlast[sQuery][iRank]["SubjectStart"],
					dBlast[sQuery][iRank]["SubjectEnd"],dBlast[sQuery][iRank]["Evalue"],
					dBlast[sQuery][iRank]["BitScore"],dContent[sQuery]
					]
				else:
					sReadQuantity=dContigs[UNASSIGNED_READS]
					tLine=[sRank,sQuery,sGlobalSample,sReadQuantity,str(iQuerySize),
					DEFAULT,DEFAULT,
					DEFAULT,DEFAULT,
					DEFAULT,sSubjectId,
					sOrganism,sSuperkingdom,
					sLineage,sDefinition,str(fFragment),
					dBlast[sQuery][iRank]["Identity"],
					str(fCover),dBlast[sQuery][iRank]["Length"],dBlast[sQuery][iRank]["Mismatch"],
					dBlast[sQuery][iRank]["GapOpen"],dBlast[sQuery][iRank]["QueryStart"],
					dBlast[sQuery][iRank]["QueryEnd"],dBlast[sQuery][iRank]["SubjectStart"],
					dBlast[sQuery][iRank]["SubjectEnd"],dBlast[sQuery][iRank]["Evalue"],
					dBlast[sQuery][iRank]["BitScore"],dContent[sQuery]
					]
				FILE.write("\t".join(tLine)+"\n")
						
# HEADER_LIST=["Hit rank","Query Seq-Id","Sample","Read quantity","Sequence length","Location","Date","Host",
# "Individual","Weight(mg)","Subject Seq-Id","Organism","SuperKingdom","Taxonomy","Hit definition","% Fragment","Identity",
# "Query cover","Alignment length","Mismatches","Gap opening","Start alignment query","End alignment query",
# "Start alignment subject","End alignment subject","E-value","Bit score","Query sequences"]

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
	dData=LoadData(sDataFile)
	dMetadata=LoadMetadata(dData)
	dLength=LoadLength(sLengthFile)
	FILE=open(BLAST_OUTPUT,"w")
	FILE.write(HEADER)
	for iIndex in range(1,iJobs+1):
		print("Working on index "+str(iIndex))
		dQuery2Content=LoadQuery(BLAST_FOLDER+"/"+BLAST_INPUT.replace(REPLACEME,str(iIndex)))
		dContigs2Sample2Quantity=LoadContigsAndQuantity(REVERSE_ASSEMBLY,dQuery2Content)
		dTaxo=LoadTaxo(BLAST_FOLDER+"/"+TAXO_FILE.replace(REPLACEME,str(iIndex)))
		dBlast=LoadBlast(BLAST_FOLDER+"/"+BLAST_FILE.replace(REPLACEME,str(iIndex)))
		WriteData(FILE,dBlast,dTaxo,dContigs2Sample2Quantity,dMetadata,dQuery2Content,dLength)
	FILE.close()
	
########################################################################    
iTime2=time.time()
iDeltaTime=iTime2-iTime1
print("Script done: "+str(iDeltaTime))
