# coding: utf-8
"""Python3.6"""
# compatibility: python2.7, python2.6

import time
from optparse import OptionParser
from CreateTable import *

sCurrentVersionScript="v1"
iTime1=time.time()
########################################################################
'''
V2-2019/11/15
Import function from basic script

V1-2019/11/06
Create tsv bilan from BlastX and BlastN data

python CreateTable.py -i JOBS -p PID
JOBS: jobs number
PID: Plaque Id
'''
########################################################################
#CONSTANT
HEADER_LIST=["Hit rank","Query Seq-Id","Sample","Read quantity","Sequence length","Location","Date","Host",
"Individual","Weight(mg)","Subject Seq-Id","Validation","Organism","SuperKingdom","Taxonomy","Hit definition","% Fragment","Identity",
"Query cover","Alignment length","Mismatches","Gap opening","Start alignment query","End alignment query",
"Start alignment subject","End alignment subject","E-value","Bit score","Query sequences"]
HEADER="\t".join(HEADER_LIST)+"\n"

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
def FusionDict(dDict1,dDict2):
	dDict={}
	for sKey in dDict1:
		dDict[sKey]=dDict1[sKey]
	for sKey in dDict2:
		dDict[sKey]=dDict2[sKey]
	return dDict
	
def FusionBlastDict(dBlastN,dBlastX):
	dDict={}
	tKey=list(set(list(dBlastN.keys())+list(dBlastX.keys())))
	
	for sQueryId in tKey:
		dDict[sQueryId]={}
		bX=True
		bN=True
		try:
			dX=dBlastX[sQueryId]
		except KeyError:
			bX=False
		try:
			dN=dBlastN[sQueryId]
		except KeyError:
			bN=False
		
		if bX and not bN:
			#Identified only with BlastX
			for iRank in dBlastX[sQueryId]:
				dDict[sQueryId][iRank]=dBlastX[sQueryId][iRank]
				dDict[sQueryId][iRank]["Confidence"]="Single (X)"
		elif not bX and bN:
			#Identified only with BlastN
			for iRank in dBlastN[sQueryId]:
				dDict[sQueryId][iRank]=dBlastN[sQueryId][iRank]
				dDict[sQueryId][iRank]["Confidence"]="Single (N)"
		else:
			#Identified with both BlastX and BlastN
			dSubjectId2MaxBitScore={}
			dSubjectId2Confidence={}
			dSubjectId2Coord={}
			for iRank in dBlastN[sQueryId]:
				sSubjectId=dBlastN[sQueryId][iRank]["SubjectId"]
				dSubjectId2MaxBitScore[sSubjectId]=float(dBlastN[sQueryId][iRank]["BitScore"])
				dSubjectId2Confidence[sSubjectId]="Single (N)"
				dSubjectId2Coord[sSubjectId]=("N",iRank)
			for iRank in dBlastX[sQueryId]:
				sSubjectId=dBlastX[sQueryId][iRank]["SubjectId"]
				fBitScore=float(dBlastX[sQueryId][iRank]["BitScore"])
				if sSubjectId in dSubjectId2MaxBitScore: #In BlastX and BlastN
					dSubjectId2Confidence[sSubjectId]="Double"
					if fBitScore>dSubjectId2MaxBitScore[sSubjectId]:
						dSubjectId2MaxBitScore[sSubjectId]=fBitScore
						dSubjectId2Coord[sSubjectId]=("X",iRank)
				else: #Only in BlastX
					dSubjectId2MaxBitScore[sSubjectId]=fBitScore
					dSubjectId2Confidence[sSubjectId]="Single (X)"
					dSubjectId2Coord[sSubjectId]=("X",iRank)
			
			#Order bitscore
			dBitScore2SubjectId={}
			for sSubjectId in dSubjectId2MaxBitScore:
				if sSubjectId=="":
					print(sQueryId)
					print(dBlastN[sQueryId])
					print(dBlastX[sQueryId])
				try:
					dBitScore2SubjectId[dSubjectId2MaxBitScore[sSubjectId]].append(sSubjectId)
				except KeyError:
					dBitScore2SubjectId[dSubjectId2MaxBitScore[sSubjectId]]=[sSubjectId]
			#Reorder equal bitScore by Double or Single for the BestHit position
			for fBitScore in sorted(dBitScore2SubjectId,reverse=True):
				if len(dBitScore2SubjectId[fBitScore])!=1:
					tKeepDouble=[X for X in dBitScore2SubjectId[fBitScore] if dSubjectId2Confidence[X]=="Double"]
					tAll=list(tKeepDouble)+[X for X in dBitScore2SubjectId[fBitScore] if dSubjectId2Confidence[X]!="Double"]
					dBitScore2SubjectId[fBitScore]=list(tAll)
				break 
			
			#Fuuuuusion !
			for fBitScore in sorted(dBitScore2SubjectId,reverse=True):
				for sSubjectId in dBitScore2SubjectId[fBitScore]:
					sConfidence=dSubjectId2Confidence[sSubjectId]
					dbCoord=dSubjectId2Coord[sSubjectId]
					iOldRank=dbCoord[1]
					if dbCoord[0]=="X":
						iNewRank=len(dDict[sQueryId])+1
						dDict[sQueryId][iNewRank]=dBlastX[sQueryId][iOldRank]
						dDict[sQueryId][iNewRank]["Confidence"]=dSubjectId2Confidence[sSubjectId]
					else:
						iNewRank=len(dDict[sQueryId])+1
						dDict[sQueryId][iNewRank]=dBlastN[sQueryId][iOldRank]
						dDict[sQueryId][iNewRank]["Confidence"]=dSubjectId2Confidence[sSubjectId]
				
			# print(dSubjectId2MaxBitScore)
			# print(dSubjectId2Confidence)
			# print(dSubjectId2Coord)
			# print(dDict[sQueryId])
					
	return dDict				

def WriteData(FILE,dBlast,dTaxo,dContigs,dMetadata,dContent,dLength):
	for sQuery in dBlast:
		for iRank in dBlast[sQuery]:
			
			if iRank==1:
				sRank=BEST_HIT
			else:
				sRank=HIT
			
			iQuerySize=len(dContent[sQuery])
			iCoverSize=int(dBlast[sQuery][iRank]["QueryEnd"])-int(dBlast[sQuery][iRank]["QueryStart"])
			fCover=round(float(iCoverSize)/iQuerySize*100,2)
			
			try:
				tSample=dContigs[sQuery]
			except KeyError:
				tSample=[sQuery.split("_")[-1]]
			
			if CONTIG in sQuery:
				sReadQuantity=sQuery.split("(")[-1].split(")")[0]
			else:
				sReadQuantity="1"
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
				dMetadata[sSample]["Weight"],sSubjectId,dBlast[sQuery][iRank]["Confidence"],
				dTaxo[sSubjectId]["Organism"],dTaxo[sSubjectId]["Superkingdom"],
				dTaxo[sSubjectId]["Lineage"],dTaxo[sSubjectId]["Definition"],str(fFragment),
				dBlast[sQuery][iRank]["Identity"],
				str(fCover),dBlast[sQuery][iRank]["Length"],dBlast[sQuery][iRank]["Mismatch"],
				dBlast[sQuery][iRank]["GapOpen"],dBlast[sQuery][iRank]["QueryStart"],
				dBlast[sQuery][iRank]["QueryEnd"],dBlast[sQuery][iRank]["SubjectStart"],
				dBlast[sQuery][iRank]["SubjectEnd"],dBlast[sQuery][iRank]["Evalue"],
				dBlast[sQuery][iRank]["BitScore"],dContent[sQuery]
				]
				FILE.write("\t".join(tLine)+"\n")

########################################################################
#MAIN
if __name__ == "__main__":
	dMetadata=LoadMetadata(sMeta)
	dLength=LoadLength(sLengthFile)
	FILE=open(BLASTALL_OUTPUT,"w")
	FILE.write(HEADER)
	for iIndex in range(1,iJobs+1):
		dQueryN2Content=LoadQuery(BLASTN_FOLDER+"/"+BLAST_INPUT.replace(REPLACEME,str(iIndex)))
		dQueryX2Content=LoadQuery(BLASTX_FOLDER+"/"+BLAST_INPUT.replace(REPLACEME,str(iIndex)))
		dQuery2Content=FusionDict(dQueryN2Content,dQueryX2Content)
		dContigs2Sample=LoadContigs(SHORTSPADES,dQuery2Content)
		dContigs2Sample=LoadContigs(SHORTFLASH,dQuery2Content,dContigs2Sample)
		dTaxoN=LoadTaxo(BLASTN_FOLDER+"/"+TAXON_FILE.replace(REPLACEME,str(iIndex)))
		dTaxoX=LoadTaxo(BLASTX_FOLDER+"/"+TAXOX_FILE.replace(REPLACEME,str(iIndex)))
		dBlastN=LoadBlast(BLASTN_FOLDER+"/"+BLASTN_FILE.replace(REPLACEME,str(iIndex)))
		dBlastX=LoadBlast(BLASTX_FOLDER+"/"+BLASTX_FILE.replace(REPLACEME,str(iIndex)))
		dTaxo=FusionDict(dTaxoN,dTaxoX)
		dBlast=FusionBlastDict(dBlastN,dBlastX)
		WriteData(FILE,dBlast,dTaxo,dContigs2Sample,dMetadata,dQuery2Content,dLength)
	FILE.close()
	
########################################################################    
iTime2=time.time()
iDeltaTime=iTime2-iTime1
print("Script done: "+str(iDeltaTime))

