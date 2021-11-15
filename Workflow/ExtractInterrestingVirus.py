# coding: utf-8
"""Python3.6"""

import time
from optparse import OptionParser
import os

sCurrentVersionScript="v1"
iTime1=time.time()
########################################################################
'''
V1-2021/10/08
Analyse bigTable to extract most interresting viruses

python ExtractInterrestingVirus.py -i INPUT -o OUTPUT -f STAT_FOLDER [-s SIZE_THRESHOLD -1 TRACK1_THRESHOLD -2 TRACK2_THRESHOLD -n NEIGHBOURS_VALUE]

INPUT: NearVANA results .tab
OUTPUT: NearVANA results .tab restricted to interresting viruses
STAT_FOLDER: The stat folder generated (likely ${PID}_Stat_BlastD/ or ${PID}_Stat_BlastN/)
SIZE_THRESHOLD: Min % of species reference genome to keep a contigs
TRACK1_THRESHOLD: Min % of reads from a sample into a contigs to keep this sample associated to the contigs
TRACK2_THRESHOLD: Min % of reads from all contigs from the same species into a sample to keep a sample associated to a contigs
NEIGHBOURS_VALUE: FLoat. During TRACK1 filtering, keep sample under TRACK1_THRESHOLD if there float*% is higher than last keeped sample
'''
########################################################################
#CONSTANT
EMPTY=""
ENDLINE="\n"
TABULATION="\t"
ZERO=0
NA="N/A"
BEST_HIT="Best hit"

STATBYCONTIGS="StatByContigs"
SAMPLE2CONTIGS="Sample2Contigs"
TAG_TSV=".tsv"

INTEREST_THRESHOLD=75
TRACK1_KEEPING_THRESHOLD=10
TRACK2_KEEPING_THRESHOLD=20
NEIGHBOURS_VALUE=1.1
UNASSIGNED_READS="UnassignedReads"

CONTENT="Content"
NAME="Name"
DEF="Def"
HOST="Host"
CHECK1="Check1"
CHECK2="Check2"

#NEARVANA TSV
COL_HITRANK=0
COL_QUERYID=1
COL_SAMPLEID=2
COL_READNUMBER=3
COL_SEQSIZE=4
COL_LOCATION=5
COL_DATE=6
COL_HOST=7
COL_HOSTNUMBER=8
COL_HOSTWEIGTH=9
COL_REFSEQID=10
COL_REFORGANISM=11
COL_REFSUPKINGDOM=12
COL_REFTAXONOMY=13
COL_REFDEFINITION=14
COL_PERCENTFRAG=15
COL_IDENTITY=16
COL_QUERYCOVER=17
COL_ALIGNMENTSIZE=18
COL_MISMATCH=19
COL_OPENGAP=20
COL_QUERYSTART=21
COL_QUERYEND=22
COL_REFSTART=23
COL_REFEND=24
COL_EXPECTATION=25
COL_BITSCORE=26
COL_QUERYSEQ=27

########################################################################
#Options
parser = OptionParser(conflict_handler="resolve")
parser.add_option("-i","--input", dest="input")
parser.add_option("-f","--folder", dest="folder")
parser.add_option("-o","--output", dest="output")
parser.add_option("-s","--size_threshold", dest="size_threshold")
parser.add_option("-1","--track1_threshold", dest="track1_threshold")
parser.add_option("-2","--track2_threshold", dest="track2_threshold")
parser.add_option("-n","--neighbours_value", dest="neighbours_value")

(options, args) = parser.parse_args()

sInput=options.input
if not sInput:
	exit("Error : no input -i defined, process broken")

sFolder=options.folder
if not sFolder:
	exit("Error : no folder -f defined, process broken")
    
sOutput=options.output
if not sOutput:
	exit("Error : no output -o defined, process broken")

sSize_Threshold=options.size_threshold
if sSize_Threshold:
    try:
        iSize_Threshold=int(sSize_Threshold)
        INTEREST_THRESHOLD=iSize_Threshold
    except KeyError:
        exit("Error : no size_threshold -s must be an integer, process broken")
else:
    print("Warning : no size_threshold -s defined, default value: {} (%)".format(INTEREST_THRESHOLD))

sTrack1_Threshold=options.track1_threshold
if sTrack1_Threshold:
    try:
        iTrack1_Threshold=int(sTrack1_Threshold)
        TRACK1_KEEPING_THRESHOLD=iTrack1_Threshold
    except KeyError:
        exit("Error : no track1_threshold -1 must be an integer, process broken")
else:
    print("Warning : no track1_threshold -1 defined, default value: {} (%)".format(TRACK1_KEEPING_THRESHOLD))

sTrack2_Threshold=options.track2_threshold
if sTrack2_Threshold:
    try:
        iTrack2_Threshold=int(sTrack2_Threshold)
        TRACK2_KEEPING_THRESHOLD=iTrack2_Threshold
    except KeyError:
        exit("Error : no track2_threshold -2 must be an integer, process broken")
else:
    print("Warning : no track21_threshold -2 defined, default value: {} (%)".format(TRACK2_KEEPING_THRESHOLD))

sNeighbours_value=options.neighbours_value
if sNeighbours_value:
    try:
        fNeighbours_value=float(sNeighbours_value)
        NEIGHBOURS_VALUE=fNeighbours_value
    except KeyError:
        exit("Error : no neighbours_value -n must be a float, process broken")
else:
    print("Warning : no neighbours_value -n defined, default value: {} (%)".format(NEIGHBOURS_VALUE))

########################################################################
#Function
def GetRelevantContigs(dDict,sPath):
    for sBestHit in dDict:
        print(sBestHit)
        tWeakContigs=[]
        for sContigId in dDict[sBestHit]:
            print("\t",sContigId,len(dDict[sBestHit][sContigId][CONTENT]))
            #Check01: Up to 10% of reads from the contigs shared by SampleId. Fuzzy limits.
            sFile="{}/{}/{}{}".format(sPath,STATBYCONTIGS,sContigId,TAG_TSV)
            
            dPercentReads2SampleId={}
            tPercent=[]
            bHeader=True   
            fHack=0.00001         
            try:
                for sNewLine in open(sFile):
                    if bHeader:
                        bHeader=False
                        continue
                    sLine=sNewLine.strip()
                    # print(sLine)
                    tLine=sLine.split(TABULATION)
                    sSampleId=tLine[0]
                    if sSampleId==UNASSIGNED_READS:
                        continue
                    fPercent=float(tLine[2])
                    try:
                        oCrash=dPercentReads2SampleId[fPercent]
                        fPercent+=fHack
                        fHack+=0.00001
                    except KeyError:
                        pass
                    dPercentReads2SampleId[fPercent]=sSampleId
                    tPercent.append(fPercent)
            except FileNotFoundError:
                print("\t /!\ Weak contigs ! Removed...")
                tWeakContigs.append(sContigId)
                continue
                
            #DEBUG
            if len(tPercent)!=len(dPercentReads2SampleId):
                exit("ERROR 96 : Divergence!")
            #/DEBUG
            tPercent=sorted(tPercent,reverse=True)
            iIndex=0
            while True:
                if tPercent[iIndex]<TRACK1_KEEPING_THRESHOLD:
                    if tPercent[iIndex]*NEIGHBOURS_VALUE<tPercent[iIndex-1]:
                        break
                iIndex+=1
                # print(iIndex,len(tPercent))
                if iIndex>=len(tPercent):
                    break
                
            dCheck1Sample={}
            for oKey in tPercent[:iIndex+1]:
                dCheck1Sample[dPercentReads2SampleId[oKey]]=round(oKey,2)
            dTemp1Sample={}
            for oKey in tPercent[iIndex+1:]:
                dTemp1Sample[dPercentReads2SampleId[oKey]]=round(oKey,2)

            tOtherSample=[dPercentReads2SampleId[X] for X in tPercent[iIndex+1:]]
            
            #Check02: Up to 20% of reads from the sample shared the same Family?Order?Species?. Hard limit.
            dCheck2Sample={}
            for sSampleId in tOtherSample:
                sFileCheck="{}/{}/{}{}".format(sPath,SAMPLE2CONTIGS,sSampleId,TAG_TSV)
                dFileContent={}
                bHeader=True
                for sNewLine in open(sFileCheck):
                    if bHeader:
                        bHeader=False
                        continue
                    sLine=sNewLine.strip()
                    tLine=sLine.split(TABULATION)
                    sNewContigId=tLine[0]
                    fReadPercent=float(tLine[2])
                    sFamily=tLine[4].replace("(PartOf)","").replace("(Maybe)","")
                    sGenus=tLine[4].replace("(PartOf)","").replace("(Maybe)","")
                    sSpecies=tLine[4].replace("(PartOf)","").replace("(Maybe)","")
                    dFileContent[sNewContigId]={"Percent":fReadPercent,"Family":sFamily,"Genus":sGenus,"Species":sSpecies}
                    
                fSumPercent=dFileContent[sContigId]["Percent"]
                sTargetKeyword=dFileContent[sContigId]["Species"]
                # print(sSampleId,sTargetKeyword,fSumPercent)
                
                for sNewContigId in dFileContent:
                    if sNewContigId==sContigId:
                        continue
                    if dFileContent[sNewContigId]["Species"]==sTargetKeyword:
                        fSumPercent+=dFileContent[sNewContigId]["Percent"]
                        # print(sNewContigId,fSumPercent)
                if fSumPercent>=TRACK2_KEEPING_THRESHOLD:
                    # tCheck2Sample.append(sSampleId)
                    dCheck2Sample[sSampleId]=fSumPercent
            
            #Update Check01
            dTemp2Sample={}
            for sSampleId in dCheck1Sample:
                sFileCheck="{}/{}/{}{}".format(sPath,SAMPLE2CONTIGS,sSampleId,TAG_TSV)
                dFileContent={}
                bHeader=True
                for sNewLine in open(sFileCheck):
                    if bHeader:
                        bHeader=False
                        continue
                    sLine=sNewLine.strip()
                    tLine=sLine.split(TABULATION)
                    sNewContigId=tLine[0]
                    fReadPercent=float(tLine[2])
                    sFamily=tLine[4].replace("(PartOf)","").replace("(Maybe)","")
                    sGenus=tLine[4].replace("(PartOf)","").replace("(Maybe)","")
                    sSpecies=tLine[4].replace("(PartOf)","").replace("(Maybe)","")
                    dFileContent[sNewContigId]={"Percent":fReadPercent,"Family":sFamily,"Genus":sGenus,"Species":sSpecies}
                    
                fSumPercent=dFileContent[sContigId]["Percent"]
                sTargetKeyword=dFileContent[sContigId]["Species"]
                # print(sSampleId,sTargetKeyword,fSumPercent)
                
                for sNewContigId in dFileContent:
                    if sNewContigId==sContigId:
                        continue
                    if dFileContent[sNewContigId]["Species"]==sTargetKeyword:
                        fSumPercent+=dFileContent[sNewContigId]["Percent"]
                        # print(sNewContigId,fSumPercent)
                dTemp2Sample[sSampleId]=fSumPercent
            
            tCheckedSample=list(dCheck1Sample.keys())+list(dCheck2Sample.keys())
            
            tRejectedSample=[]
            for sSampleId in dDict[sBestHit][sContigId][CONTENT]:
                if sSampleId not in tCheckedSample:
                    tRejectedSample.append(sSampleId)
            for sSampleId in tRejectedSample:
                del dDict[sBestHit][sContigId][CONTENT][sSampleId]
                
            for sSampleId in dCheck1Sample:
                dDict[sBestHit][sContigId][CONTENT][sSampleId][CHECK1]=dCheck1Sample[sSampleId]
                dDict[sBestHit][sContigId][CONTENT][sSampleId][CHECK2]=dTemp2Sample[sSampleId]
            for sSampleId in dCheck2Sample:
                dDict[sBestHit][sContigId][CONTENT][sSampleId][CHECK1]=dTemp1Sample[sSampleId]
                dDict[sBestHit][sContigId][CONTENT][sSampleId][CHECK2]=dCheck2Sample[sSampleId]
            
            print("\t '->",sContigId,len(dDict[sBestHit][sContigId][CONTENT]),"->",len(dCheck1Sample),"+",len(dCheck2Sample))
            # print(dDict[sBestHit][sContigId])
        
        for sContigId in tWeakContigs:
            del dDict[sBestHit][sContigId]
    
    return dDict
            
def GetInterestingContigs(sPath):
    dDict={}
    bHeader=True
    for sNewLine in open(sPath):
        if bHeader:
            bHeader=False
            continue
        sLine=sNewLine.strip()
        tLine=sLine.split(TABULATION)
        sBestHit=tLine[COL_HITRANK]
        if sBestHit!=BEST_HIT:
            continue
        try:
            oCrash=dDict[sBestHit]
        except KeyError:
            dDict[sBestHit]={}
        oFrag=tLine[COL_PERCENTFRAG]
        sContigId=tLine[COL_QUERYID]
        if oFrag==NA:
            continue
        sSampleId=tLine[COL_SAMPLEID]
        if sSampleId==UNASSIGNED_READS:
            continue
        fFrag=float(tLine[COL_PERCENTFRAG])
        if fFrag>=INTEREST_THRESHOLD:
            sHost=tLine[COL_HOST]
            sName=tLine[COL_REFORGANISM]
            sDef=tLine[COL_REFDEFINITION]
            # print(sBestHit,sContigId,sSampleId,sHost,sName,sDef,fFrag)     
            try:
                oCrash=dDict[sBestHit][sContigId]
                dDict[sBestHit][sContigId][CONTENT][sSampleId]={HOST:sHost,CHECK1:0.0,CHECK2:0.0}
            except KeyError:
                dDict[sBestHit][sContigId]={NAME:sName,DEF:sDef,CONTENT:{sSampleId:{HOST:sHost,CHECK1:0.0,CHECK2:0.0}}}
    return dDict
        
def RewriteTable(dDict,sPathRef,sPathNew):
    with open(sPathNew,"w") as oOutputFile:
        sHeader="Plate\tQuery Seq-Id\tSample\tRead quantity\tSequence length\t%Reads/contigs\t%Reads-species/sample\tLocation\tDate\tHost\tIndividual\tWeight(mg)\tSubject Seq-Id\tOrganism\tSuperKingdom\tTaxonomy\tHit definition\t% Fragment\tIdentity\tQuery cover\tAlignment length\tMismatches\tGap opening\tStart alignment query\tEnd alignment queryStart alignment subject\tEnd alignment subject\tE-value\tBit score\tQuery sequences\n"
        oOutputFile.write(sHeader)
        with open(sPathRef) as oInputFile:
            bHeader=True
            for sNewLine in oInputFile:
                if bHeader:
                    bHeader=False
                    continue
                tLine=sNewLine.split(TABULATION)
                sBestHit=tLine[COL_HITRANK]
                if sBestHit!=BEST_HIT:
                    continue
                sContigId=tLine[COL_QUERYID]
                sSampleId=tLine[COL_SAMPLEID]
                try:
                    oCrash=dDict[sBestHit][sContigId][CONTENT][sSampleId]
                    # print("Relevant line: {} {} {}".format(sBestHit,sContigId,sSampleId))
                except KeyError:
                    continue
                oCheck1=dDict[sBestHit][sContigId][CONTENT][sSampleId][CHECK1]
                oCheck2=dDict[sBestHit][sContigId][CONTENT][sSampleId][CHECK2]
                sCorrectedLine=TABULATION.join(tLine[:5])+"\t{}\t{}\t".format(oCheck1,oCheck2)+TABULATION.join(tLine[5:])
                oOutputFile.write(sCorrectedLine)

########################################################################
#MAIN
if __name__ == "__main__":
    dTarget=GetInterestingContigs(sInput)
    dTarget=GetRelevantContigs(dTarget,sFolder)
    RewriteTable(dTarget,sInput,sOutput)
    
    
########################################################################    
iTime2=time.time()
iDeltaTime=iTime2-iTime1
print("Script done: "+str(iDeltaTime))

