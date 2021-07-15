# coding: utf-8
"""Python3.6"""

import time
from optparse import OptionParser
import os

sCurrentVersionScript="v19"
iTime1=time.time()
########################################################################
'''
V19-2021/03/03
Adaptation to RNA_Seq: No Sample!

V18-2021/02/02
Retrieve Sample from Metadata in order to show empty sample
V17-2021/02/01
Filter FILTER_MINREAD2ASSIGN to reject contigs assignation to sample if reads from that sample are equal or less
that the filter.
V16-2021/01/14
Hard exception to FILTER_READ_REPRESENTATIVITE : at least 1000 reads is a confident clue of a viruses
V15-2021/01/12
Change FILTER_READ_REPRESENTATIVITE to match a ratio. Default 1.0
A contig is associated to a sample only if reads from that sample compared to all reads from the sample exceed this ratio.
V14-2020/09/25
Add 3 columns (Family,Genus,Species) into Contigs2Reads
Newfolder ByContigs, 1 file/contigs, repartition of reads onto all sample, with the LIMITREADS effect
V13-2020/07/10
Add filter to remove Sample from contig if Sample Reads < -l reads_limit
V12-2020/07/09
Add ContigsName, ContigsSize, QueryRef for stat grouped by virus (sample scale)
V11-2020/07/03
Complete rework : get reads AND contigs information
V10-2020/03/24
Compatibility NearVANA only
Need ReverseAssembly file for Megahit as argument
V9-2020/03/23
Add: Bilan for all sample merged
Fix: Count Read only one time with NearVANA results redundancy
V8-2019/11/19
Add: Bilan-> BilanByGroup,BilanBySample
V7-2019/11/18
Change: Use the VMR file from ICTV (include genome refseq, change column position)
V6-2019/11/07
Fix: Unknown sequences no more researched because data filtered on viruses
Add Research fo date format
Add normalized data on Individuals, Weigth and Weigth/Individual
V5-2019/10/16
Add filter based on Read size.
V4-2019/10/02
Too much "None" assignation -> If ICTV fail, use Blast taxo. If Blast taxo fail, unclassified
V3-2019/09/27
Don't take into account data without "Viruses" Kingdom.
V2-2019/09/26
Better folder ramification to handle new stat file
V1-2019/09/10

python Extract4Stat.py -i INPUT -o OUTPUT -v ICTV -r REVERSE [-s SIZETHRESHOLD -q QUANTITYTHRESOLD -h HARDLIMIT
                        -t REJECTED -n NOISE]
INPUT: NearVANA tsv table
OUTPUT: output folder results
ICTV: ictv taxonomy tsv table
REVERSE: ReverseAssembly file
SIZETHRESHOLD: Minimum size to take a contig into account
QUANTITYTHRESOLD: Minimum quantity of reads to take a contig into account
HARDLIMIT: Minimum value to counter QUANTITYTHRESHOLD
REJECTED: List of taxo rejected, separated by "-". Space " " must be replace by underscore "_"
NOISE: Part of Contigs below or equal to this value are ignored 
'''
########################################################################
#CONSTANT
BEST_HIT="Best hit"
TARGET="Viruses"

FILTER_CONTIG_SIZE=0
FILTER_READ_QUANTITY=300
FILTER_READ_HARDLIMIT=1000
FILTER_NOISE=0

INTEGER=["0","1","2","3","4","5","6","7","8","9"]
SPACE=" "
UNDERSCORE="_"

REJECTED_SEPARATOR="-"
REJECTED_TAXO=[]

#NEARVANA RNAseq TSV
COL_HITRANK=0
COL_QUERYID=1
COL_READNUMBER=2
COL_SEQSIZE=3
COL_REFSEQID=4
COL_REFORGANISM=5
COL_REFSUPKINGDOM=6
COL_REFTAXONOMY=7
COL_REFDEFINITION=8
COL_PERCENTFRAG=9
COL_IDENTITY=10
COL_QUERYCOVER=11
COL_ALIGNMENTSIZE=12
COL_MISMATCH=13
COL_OPENGAP=14
COL_QUERYSTART=15
COL_QUERYEND=16
COL_REFSTART=17
COL_REFEND=18
COL_EXPECTATION=19
COL_BITSCORE=20
COL_QUERYSEQ=21

UNCLASSIFIED="unclassified"
UNKNOWN="unknown"
UNASSIGNED="unassigned"

FAMILY="Family"
GENUS="Genus"
SPECIES="Species"

READS="Reads"
CONTIGS="Contigs"
TAXO="Taxo"

SIZE="Size"
REF="Ref"

#ICTV TSV
TAXONOMY=["Realm","Subrealm","Kingdom","Subkingdom","Phylum","Subphylum","Class","Subclass","Order","Suborder","Family","Subfamily","Genus","Subgenus","Species"]
KEYWORD_FAMILY="Family"
KEYWORD_GENUS="Genus"
KEYWORD_SPECIES="Species"

ABUNDANCY_THRESHOLD=10

MAYBE="(Maybe)"
FAMILY_THRESHOLD=00
GENUS_THRESHOLD=30
SPECIES_THRESHOLD=80
PARTOF="(Partof)"
COVER_THRESHOLD=66

FOLDER_FAMILYGLOBAL="StatByFamily"
FOLDER_FAMILY="ByFamily"
FOLDER_GENUS="ByGenus"
FOLDER_SPECIES="BySpecies"
FOLDER_PLATE="StatPlate"

LOG="options.log"

MASTER_SAMPLE_TAG="------ Get Sample list ------"

########################################################################
#Options
parser = OptionParser(conflict_handler="resolve")
parser.add_option("-i","--input", dest="input")
parser.add_option("-o","--output", dest="output")
parser.add_option("-v","--taxo_ictv", dest="taxo_ictv")
parser.add_option("-q","--quantity_treshold", dest="quantity_threshold")
parser.add_option("-s","--size_threshold", dest="size_threshold")
parser.add_option("-t","--rejected_taxo", dest="rejected_taxo")
parser.add_option("-r","--reverse_assembly", dest="reverse_assembly")
parser.add_option("-n","--noise", dest="noise")
parser.add_option("-h","--hardlimit", dest="hardlimit")

(options, args) = parser.parse_args()

sInput=options.input
if not sInput:
	exit("Error : no input -i defined, process broken")

sOutput=options.output
if not sOutput:
	exit("Error : no output folder -o defined, process broken")

sTaxo=options.taxo_ictv
if not sTaxo:
	exit("Error : no taxo_ictv -v defined, process broken")
    
sReverseAssembly=options.reverse_assembly
if not sReverseAssembly:
    exit("Error : no reverse_assembly -r defined, process broken")    

sQuantityThreshold=options.quantity_threshold
if not sQuantityThreshold:
    print("Warning : no quantity_threshold -q defined, default value : "+str(FILTER_READ_QUANTITY))
    sQuantityThreshold=str(FILTER_READ_QUANTITY)
try:
    iQuantityThreshold=int(sQuantityThreshold)
    FILTER_READ_QUANTITY=iQuantityThreshold
except ValueError:
    exit("Error : no quantity_threshold -q must be an integer, process broken")

sSizeThreshold=options.size_threshold
if not sSizeThreshold:
    print("Warning : no size_threshold -s defined, default value : "+str(FILTER_CONTIG_SIZE))
    sSizeThreshold=str(FILTER_CONTIG_SIZE)
try:
    iSizeThreshold=int(sSizeThreshold)
    FILTER_CONTIG_SIZE=iSizeThreshold
except ValueError:
    exit("Error : no size_threshold -s must be an integer, process broken")
    
sRejected=options.rejected_taxo
if not sRejected:
    print("Warning : no rejected_taxo -t defined, all taxonomy will be considered")
    sRejected=""
else:
    sRejected=sRejected.replace(UNDERSCORE,SPACE)
    tRejected=sRejected.split(REJECTED_SEPARATOR)
    REJECTED_TAXO=tRejected

sNoise=options.noise
if not sNoise:
    print("Warning : no noise -n defined, default value : "+str(FILTER_NOISE))
    sNoise=str(FILTER_NOISE)
try:
    iNoise=int(sNoise)
    FILTER_NOISE=iNoise
except ValueError:
    exit("Error : noise -n must be an integer, process broken")

sHardLimit=options.hardlimit
if not sHardLimit:
    print("Warning : no hardlimit -h defined, default value : "+str(FILTER_READ_HARDLIMIT))
    sHardLimit=str(FILTER_READ_HARDLIMIT)
try:
    iHardLimit=int(sHardLimit)
    FILTER_READ_HARDLIMIT=iHardLimit
except ValueError:
    exit("Error : hardlimit -h must be an integer, process broken")

########################################################################
#Function
def StoreLog():
    FILE=open(sOutput+"/"+LOG,"w")
    FILE.write("script_version: "+sCurrentVersionScript+"\n")
    FILE.write("input: "+sInput+"\n")
    FILE.write("output: "+sOutput+"\n")
    FILE.write("taxo_ictv: "+sTaxo+"\n")
    FILE.write("quantity_treshold: "+sQuantityThreshold+"\n")
    FILE.write("hardlimit: "+sHardLimit+"\n")
    FILE.write("size_threshold: "+sSizeThreshold+"\n")
    FILE.write("rejected_taxo: "+sRejected+"\n")
    FILE.write("reverse_assembly: "+sReverseAssembly+"\n")
    FILE.write("reads_noise_threshold: "+sNoise+"\n")
    FILE.write("PartOf_coverThresold: "+str(COVER_THRESHOLD)+"\n")
    FILE.write("Maybe_familyThreshold: "+str(FAMILY_THRESHOLD)+"\n")
    FILE.write("Maybe_genusThreshold: "+str(GENUS_THRESHOLD)+"\n")
    FILE.write("Maybe_speciesThreshold: "+str(SPECIES_THRESHOLD)+"\n")
    FILE.close()

def CreateFolder(sPath):
    tFolder=[sPath,"{}/{}".format(sPath,FOLDER_FAMILYGLOBAL),
            "{}/{}".format(sPath,FOLDER_PLATE),
            ]
    for sFolder in tFolder:  
        try:
            os.mkdir(sFolder)
        except OSError or FileExistsError:
            print("Warning: Can not create output folder {}".format(sFolder))

def LoadICTV(sPath):
    bHeader=True
    dResult={}
    dKeyword2ColId={}
    iLine=0
    for sNewLine in open(sPath):
        iLine+=1
        sLine=sNewLine.strip()
        tLine=sLine.split("\t")
        if bHeader:
            for iColIndex in range(len(tLine)):
                if tLine[iColIndex] in TAXONOMY:
                    iTaxIndex=TAXONOMY.index(tLine[iColIndex])
                    dKeyword2ColId[TAXONOMY[iTaxIndex]]=iColIndex
            bHeader=False
            continue
        sFamily=tLine[dKeyword2ColId[KEYWORD_FAMILY]]
        if sFamily!="" and sFamily!=" ":
            dResult[sFamily]=KEYWORD_FAMILY
        sGenus=tLine[dKeyword2ColId[KEYWORD_GENUS]]
        if sGenus!="" and sGenus!=" ":
            dResult[sGenus]=KEYWORD_GENUS
        sSpecies=tLine[dKeyword2ColId[KEYWORD_SPECIES]]
        if sSpecies!="" and sSpecies!=" ":
            dResult[sSpecies]=KEYWORD_SPECIES
    return dResult
    
def LoadData(sPath):
    dResult={}
    bHeader=True
    for sNewLine in open(sPath):
        if bHeader:
            bHeader=False
            continue
        sLine=sNewLine.strip()
        tLine=sLine.split("\t")
        if tLine[COL_HITRANK]==BEST_HIT and tLine[COL_REFSUPKINGDOM]==TARGET \
            and int(tLine[COL_SEQSIZE])>=FILTER_CONTIG_SIZE and int(tLine[COL_READNUMBER])>=FILTER_READ_QUANTITY:
            sTaxo=tLine[COL_REFTAXONOMY]
            bBreak=False
            for sRejected in REJECTED_TAXO:
                if sRejected in sTaxo:
                    bBreak=True
                    break
            if bBreak:
                continue
            sUniqSeqName=tLine[COL_QUERYID]
            dbLine=tuple(tLine)
            dResult[sUniqSeqName]=dbLine
    return dResult    

def ReverseTableDict(dDict):    
    dResult={}
    for sKey in dDict:
        for sValue in dDict[sKey]:
            dResult[sValue]=sKey
    return dResult

def LoadReverseAssembly(dData,sFile):
    dResult={}
    for sKey in dData:
        sContigName=sKey
        try:
            oCrash=dResult[sContigName]
        except KeyError:
            dResult[sContigName]=0
   
    for sNewLine in open(sFile):
        tLine=sNewLine.split("\t")
        sContigName=tLine[1]
        try:
            oCrash=dResult[sContigName]
        except KeyError:
            continue
        sRead=tLine[0]
        dResult[sContigName]+=1

    if FILTER_NOISE!=0:
        tRemoveContig=[]
        for sContigName in dResult:
            #Remove reads quantity below T- sample
            if dResult[sContigName]<=FILTER_NOISE:
                tRemoveContig.append(sContigName)
        for sContigName in tRemoveContig:
            del dResult[sContigName]
            
    if FILTER_READ_HARDLIMIT!=0:
        tRemoveContig=[]
        for sContigName in dResult:
            #Remove reads quantity below T- sample
            if dResult[sContigName]<=FILTER_READ_HARDLIMIT:
                tRemoveContig.append(sContigName)
        for sContigName in tRemoveContig:
            del dResult[sContigName]
            
    return dResult

def GetTaxo(dData):
    dResult={}
    for sItem in dData:        
        dbLine=dData[sItem]
        sContigName=dbLine[COL_QUERYID]
        sTaxonomy=dbLine[COL_REFTAXONOMY]
        fIdentity=float(dbLine[COL_IDENTITY])
        fCover=float(dbLine[COL_QUERYCOVER])
        
        sTaxonomy=sTaxonomy.replace("; ",";")
        tTaxonomy=sTaxonomy.split(";")
        while len(tTaxonomy[-1])==0:
            tTaxonomy=tTaxonomy[:-1]
        sFamily=None
        sGenus=None
        sSpecies=None
        for sKey in tTaxonomy:
            try:
                sTaxo=dICTV[sKey]
                if sTaxo==KEYWORD_FAMILY:
                    sFamily=sKey
                elif sTaxo==KEYWORD_GENUS:
                    sGenus=sKey
                elif sTaxo==KEYWORD_SPECIES:
                    sSpecies=sKey
            except KeyError:
                continue
                
        if sFamily is None and sSpecies is None and sGenus is None:
            sFamily=sGenus=sSpecies=tTaxonomy[-1]
        elif sFamily is None and (sGenus is not None or sSpecies is not None):
            sFamily="{} {}".format(UNASSIGNED,FAMILY)
        
        if sGenus is None and sSpecies is None:
            if UNCLASSIFIED in tTaxonomy[-1]:
                sGenus=tTaxonomy[-1]
                sSpecies=tTaxonomy[-1]
            else:
                sGenus="{} {}".format(UNKNOWN,sFamily)
                sSpecies="{} {}".format(UNKNOWN,sFamily)
        elif sGenus is None and sSpecies is not None:
            if UNCLASSIFIED in tTaxonomy[-1]:
                sGenus=tTaxonomy[-1]
            else:
                sGenus="{} {}".format(UNKNOWN,sFamily)
        
        if sSpecies is None:
            if UNCLASSIFIED in tTaxonomy[-1]:
                sSpecies=tTaxonomy[-1]
            else:
                sSpecies="{} {}".format(UNKNOWN,sGenus)
        
        if sFamily is None or sGenus is None or sSpecies is None:
            print(tTaxonomy)
            exit()
        
        if abs(fCover)<COVER_THRESHOLD:
            sFamily=PARTOF+sFamily
            sGenus=PARTOF+sGenus
            sSpecies=PARTOF+sSpecies
        if fIdentity<FAMILY_THRESHOLD:
            sFamily=MAYBE+sFamily
        if fIdentity<GENUS_THRESHOLD:
            sGenus=MAYBE+sGenus
        if fIdentity<SPECIES_THRESHOLD:
            sSpecies=MAYBE+sSpecies
        
        dResult[sContigName]={FAMILY:sFamily,GENUS:sGenus,SPECIES:sSpecies}
    
    return dResult

def ProcessStatAllSample(sOutput,dSeq2Quantity,dContig2Taxo):
    dTaxoFamily={}
    dTaxoGenus={}
    dTaxoSpecies={}
    iSumContigs=0
    iSumReads=0
    for sSeq in dSeq2Quantity:
        iSumContigs+=1
        sFamily=dContig2Taxo[sSeq][FAMILY]
        sGenus=dContig2Taxo[sSeq][GENUS]
        sSpecies=dContig2Taxo[sSeq][SPECIES]
        try:
            oCrash=dTaxoFamily[sFamily]
        except KeyError:
            dTaxoFamily[sFamily]={READS:0,CONTIGS:0}
        try:
            oCrash=dTaxoGenus[sGenus]
        except KeyError:
            dTaxoGenus[sGenus]={READS:0,CONTIGS:0}
        try:
            oCrash=dTaxoSpecies[sSpecies]
        except KeyError:
            dTaxoSpecies[sSpecies]={READS:0,CONTIGS:0}
        dTaxoFamily[sFamily][CONTIGS]+=1
        dTaxoGenus[sGenus][CONTIGS]+=1
        dTaxoSpecies[sSpecies][CONTIGS]+=1        
        dTaxoFamily[sFamily][READS]+=dSeq2Quantity[sSeq]
        dTaxoGenus[sGenus][READS]+=dSeq2Quantity[sSeq]
        dTaxoSpecies[sSpecies][READS]+=dSeq2Quantity[sSeq]
        iSumReads+=dSeq2Quantity[sSeq]
    FILE_FAMILY=open(sOutput+"/"+FOLDER_PLATE+"/ByFamily.tsv","w")
    FILE_FAMILY.write("Family\tReads\tContigs\t%Reads(/All)\t%Contigs(/All)\n")
    for sFamily in sorted(dTaxoFamily):
        fReads=round(float(dTaxoFamily[sFamily][READS])/iSumReads*100,2)
        fContigs=round(float(dTaxoFamily[sFamily][CONTIGS])/iSumContigs*100,2)
        FILE_FAMILY.write("{}\t{}\t{}\t{}\t{}\n".format(sFamily,dTaxoFamily[sFamily][READS],dTaxoFamily[sFamily][CONTIGS],fReads,fContigs))
    FILE_FAMILY.close()
    FILE_GENUS=open(sOutput+"/"+FOLDER_PLATE+"/ByGenus.tsv","w")
    FILE_GENUS.write("Family\tReads\tContigs\t%Reads(/All)\t%Contigs(/All)\n")
    for sGenus in sorted(dTaxoGenus):
        fReads=round(float(dTaxoGenus[sGenus][READS])/iSumReads*100,2)
        fContigs=round(float(dTaxoGenus[sGenus][CONTIGS])/iSumContigs*100,2)
        FILE_GENUS.write("{}\t{}\t{}\t{}\t{}\n".format(sGenus,dTaxoGenus[sGenus][READS],dTaxoGenus[sGenus][CONTIGS],fReads,fContigs))
    FILE_GENUS.close()    
    FILE_SPECIES=open(sOutput+"/"+FOLDER_PLATE+"/BySpecies.tsv","w")
    FILE_SPECIES.write("Family\tReads\tContigs\t%Reads(/All)\t%Contigs(/All)\n")
    for sSpecies in sorted(dTaxoSpecies):
        fReads=round(float(dTaxoSpecies[sSpecies][READS])/iSumReads*100,2)
        fContigs=round(float(dTaxoSpecies[sSpecies][CONTIGS])/iSumContigs*100,2)
        FILE_SPECIES.write("{}\t{}\t{}\t{}\t{}\n".format(sSpecies,dTaxoSpecies[sSpecies][READS],dTaxoSpecies[sSpecies][CONTIGS],fReads,fContigs))
    FILE_SPECIES.close()        
        
def ProcessStatGroupedByVirus(sOutput,dSeq2Quantity,dContig2Taxo,dContig2Attribut):
    dFamily2Reads={}
    dFamily2Contigs={}
    iFullReads=0
    
    for sSeq in dSeq2Quantity:
        sFamily=dContig2Taxo[sSeq][FAMILY]
        sFamily=sFamily.replace(MAYBE,"").replace(PARTOF,"")
        try:
            oCrash=dFamily2Reads[sFamily]
        except KeyError:
            dFamily2Reads[sFamily]={}
        try:
            oCrash=dFamily2Contigs[sFamily]
        except KeyError:
            dFamily2Contigs[sFamily]={}
        dFamily2Reads[sFamily][sSeq]=dSeq2Quantity[sSeq]
        dFamily2Contigs[sFamily][sSeq]=True
        iFullReads+=dSeq2Quantity[sSeq]
    
    for sFamily in dFamily2Contigs:
        FILE_SAMPLE=open(sOutput+"/"+FOLDER_FAMILYGLOBAL+"/"+sFamily.replace(" ","_")+".tsv","w")
        FILE_SAMPLE.write("Contig\tReads\t%Reads(/Plates)\tContigsSize\tRefId\n")
        for sContig in sorted(dFamily2Contigs[sFamily]):
            fReads=round(float(dFamily2Reads[sFamily][sContig])/iFullReads*100,2)
            FILE_SAMPLE.write("{}\t{}\t{}\t{}\t{}\n".format(
                                        sContig,dFamily2Reads[sFamily][sContig],
                                        fReads,dContig2Attribut[sContig][SIZE],
                                        dContig2Attribut[sContig][REF]
                                        ))
        FILE_SAMPLE.close()
 
def ProcessStatByGroup(sOutput,dSeq2Sample2Quantity,dContig2Taxo,dGroup2Sample):
    dSample2Group=ReverseTableDict(dGroup2Sample)
    dTaxo2Group2Data={FAMILY:{},GENUS:{},SPECIES:{}}
    dGroup2SumContigs={}
    dGroup2SumReads={}
    for sSeq in dSeq2Sample2Quantity:
        dGroup2Bool={}
        sFamily=dContig2Taxo[sSeq][FAMILY]
        sGenus=dContig2Taxo[sSeq][GENUS]
        sSpecies=dContig2Taxo[sSeq][SPECIES]
        for sSample in dSeq2Sample2Quantity[sSeq]:
            dbGroup=dSample2Group[sSample]
            try:
                oCrash=dGroup2SumContigs[dbGroup]
            except KeyError:
                dGroup2SumContigs[dbGroup]=0
            try:
                oCrash=dGroup2SumReads[dbGroup]
            except KeyError:
                dGroup2SumReads[dbGroup]=0
            dGroup2SumContigs[dbGroup]+=1
            dGroup2SumReads[dbGroup]+=dSeq2Sample2Quantity[sSeq][sSample]
            try:
                oCrash=dTaxo2Group2Data[FAMILY][dbGroup]
            except KeyError:
                dTaxo2Group2Data[FAMILY][dbGroup]={}
            try:
                oCrash=dTaxo2Group2Data[FAMILY][dbGroup][sFamily]
            except KeyError:
                dTaxo2Group2Data[FAMILY][dbGroup][sFamily]={READS:0,CONTIGS:0}
            try:
                oCrash=dTaxo2Group2Data[GENUS][dbGroup]
            except KeyError:
                dTaxo2Group2Data[GENUS][dbGroup]={}
            try:
                oCrash=dTaxo2Group2Data[GENUS][dbGroup][sGenus]
            except KeyError:
                dTaxo2Group2Data[GENUS][dbGroup][sGenus]={READS:0,CONTIGS:0}
            try:
                oCrash=dTaxo2Group2Data[SPECIES][dbGroup]
            except KeyError:
                dTaxo2Group2Data[SPECIES][dbGroup]={}
            try:
                oCrash=dTaxo2Group2Data[SPECIES][dbGroup][sSpecies]
            except KeyError:
                dTaxo2Group2Data[SPECIES][dbGroup][sSpecies]={READS:0,CONTIGS:0}
            dGroup2Bool[dbGroup]=True
            dTaxo2Group2Data[FAMILY][dbGroup][sFamily][READS]+=dSeq2Sample2Quantity[sSeq][sSample]
            dTaxo2Group2Data[GENUS][dbGroup][sGenus][READS]+=dSeq2Sample2Quantity[sSeq][sSample]
            dTaxo2Group2Data[SPECIES][dbGroup][sSpecies][READS]+=dSeq2Sample2Quantity[sSeq][sSample]
        for dbGroup in dGroup2Bool:
            dTaxo2Group2Data[FAMILY][dbGroup][sFamily][CONTIGS]+=1
            dTaxo2Group2Data[GENUS][dbGroup][sGenus][CONTIGS]+=1
            dTaxo2Group2Data[SPECIES][dbGroup][sSpecies][CONTIGS]+=1
    
    for dbGroup in dGroup2SumContigs:
        dTarget={FAMILY:FOLDER_FAMILY,GENUS:FOLDER_GENUS,SPECIES:FOLDER_SPECIES}
        iSumContigs=dGroup2SumContigs[dbGroup]
        iSumReads=dGroup2SumReads[dbGroup]
        for sTaxo in dTarget:
            FILE_GROUP=open(sOutput+"/"+FOLDER_GROUP+"/"+dTarget[sTaxo]+"/"+FolderNameFormat(dbGroup)+".tsv","w")
            FILE_GROUP.write(sTaxo+"\tReads\tContigs\t%Reads(/Group)\t%Contigs(/Group)\n")
            for sItem in sorted(dTaxo2Group2Data[sTaxo][dbGroup]):
                iGroupContigs=dTaxo2Group2Data[sTaxo][dbGroup][sItem][CONTIGS]
                iGroupReads=dTaxo2Group2Data[sTaxo][dbGroup][sItem][READS]
                fContigs=round(float(iGroupContigs)/iSumContigs*100,2)
                fReads=round(float(iGroupReads)/iSumReads*100,2)
                FILE_GROUP.write("{}\t{}\t{}\t{}\t{}\n".format(sItem,iGroupReads,iGroupContigs,fReads,fContigs))
            FILE_GROUP.close()

def ProcessStatBySample(sOutput,dSeq2Sample2Quantity,dContig2Taxo,tSampleList):
    dTaxo2Sample2Data={FAMILY:{},GENUS:{},SPECIES:{}}
    dSample2SumContigs={}
    dSample2SumReads={}
    for sSeq in dSeq2Sample2Quantity:
        dSample2Bool={}
        sFamily=dContig2Taxo[sSeq][FAMILY]
        sGenus=dContig2Taxo[sSeq][GENUS]
        sSpecies=dContig2Taxo[sSeq][SPECIES]
        for sSample in dSeq2Sample2Quantity[sSeq]:
            try:
                oCrash=dSample2SumContigs[sSample]
            except KeyError:
                dSample2SumContigs[sSample]=0
            try:
                oCrash=dSample2SumReads[sSample]
            except KeyError:
                dSample2SumReads[sSample]=0
            dSample2SumContigs[sSample]+=1
            dSample2SumReads[sSample]+=dSeq2Sample2Quantity[sSeq][sSample]
            try:
                oCrash=dTaxo2Sample2Data[FAMILY][sSample]
            except KeyError:
                dTaxo2Sample2Data[FAMILY][sSample]={}
            try:
                oCrash=dTaxo2Sample2Data[FAMILY][sSample][sFamily]
            except KeyError:
                dTaxo2Sample2Data[FAMILY][sSample][sFamily]={READS:0,CONTIGS:0}
            try:
                oCrash=dTaxo2Sample2Data[GENUS][sSample]
            except KeyError:
                dTaxo2Sample2Data[GENUS][sSample]={}
            try:
                oCrash=dTaxo2Sample2Data[GENUS][sSample][sGenus]
            except KeyError:
                dTaxo2Sample2Data[GENUS][sSample][sGenus]={READS:0,CONTIGS:0}
            try:
                oCrash=dTaxo2Sample2Data[SPECIES][sSample]
            except KeyError:
                dTaxo2Sample2Data[SPECIES][sSample]={}
            try:
                oCrash=dTaxo2Sample2Data[SPECIES][sSample][sSpecies]
            except KeyError:
                dTaxo2Sample2Data[SPECIES][sSample][sSpecies]={READS:0,CONTIGS:0}
            dSample2Bool[sSample]=True
            dTaxo2Sample2Data[FAMILY][sSample][sFamily][READS]+=dSeq2Sample2Quantity[sSeq][sSample]
            dTaxo2Sample2Data[GENUS][sSample][sGenus][READS]+=dSeq2Sample2Quantity[sSeq][sSample]
            dTaxo2Sample2Data[SPECIES][sSample][sSpecies][READS]+=dSeq2Sample2Quantity[sSeq][sSample]
        for sSample in dSample2Bool:
            dTaxo2Sample2Data[FAMILY][sSample][sFamily][CONTIGS]+=1
            dTaxo2Sample2Data[GENUS][sSample][sGenus][CONTIGS]+=1
            dTaxo2Sample2Data[SPECIES][sSample][sSpecies][CONTIGS]+=1
    
    for sSample in tSampleList: #dSample2SumContigs:
        dTarget={FAMILY:FOLDER_FAMILY,GENUS:FOLDER_GENUS,SPECIES:FOLDER_SPECIES}
        try:
            iSumContigs=dSample2SumContigs[sSample]
            iSumReads=dSample2SumReads[sSample]
        except KeyError:
            pass
        for sTaxo in dTarget:
            FILE_SAMPLE=open(sOutput+"/"+FOLDER_SAMPLE+"/"+dTarget[sTaxo]+"/"+sSample+".tsv","w")
            FILE_SAMPLE.write(sTaxo+"\tReads\tContigs\t%Reads(/Sample)\t%Contigs(/Sample)\n")
            try:
                for sItem in sorted(dTaxo2Sample2Data[sTaxo][sSample]):
                    iSampleContigs=dTaxo2Sample2Data[sTaxo][sSample][sItem][CONTIGS]
                    iSampleReads=dTaxo2Sample2Data[sTaxo][sSample][sItem][READS]
                    fContigs=round(float(iSampleContigs)/iSumContigs*100,2)
                    fReads=round(float(iSampleReads)/iSumReads*100,2)
                    FILE_SAMPLE.write("{}\t{}\t{}\t{}\t{}\n".format(sItem,iSampleReads,iSampleContigs,fReads,fContigs))
            except KeyError:
                pass
            FILE_SAMPLE.close()
            
def GetContigAttribut(dData):
    dResult={}
    for sItem in dData:        
        dbLine=dData[sItem]
        sContigName=dbLine[COL_QUERYID]
        sSize=dbLine[COL_SEQSIZE]
        sRef=dbLine[COL_REFSEQID]
        dResult[sContigName]={SIZE:sSize,REF:sRef}
    
    return dResult

def WriteReadsByContigs(dDict,dTaxo):
    dTemp={}
    dSize={}
    iSumReads=0
    for sContigName in dDict:
        iReadsQuantity=int(sContigName.split("(")[-1].split(")")[0])
        sContigShortName=sContigName #"_".join(sContigName.split("_")[:-1])
        dSize[sContigShortName]=int(dDict[sContigName][COL_SEQSIZE])
        try:
            if sContigShortName not in dTemp[iReadsQuantity]:
                dTemp[iReadsQuantity].append(sContigShortName)
                iSumReads+=iReadsQuantity
        except KeyError:
            dTemp[iReadsQuantity]=[sContigShortName]
            iSumReads+=iReadsQuantity
    FILE=open(sOutput+"/Contigs2ReadsQuantity.tsv","w")
    FILE.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("ContigId","Contig size","Reads quantity","Reads %(/Plate)","Cumulated %","Family","Genus","Species"))
    fCumulate=0.0
    for iValue in sorted(dTemp):
        for sContigName in sorted(dTemp[iValue]):
            fCumulate+=round(float(iValue)/iSumReads*100,4)
            FILE.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(sContigName,dSize[sContigName],iValue,round(float(iValue)/iSumReads*100,4),fCumulate,dTaxo[sContigName][FAMILY],dTaxo[sContigName][GENUS],dTaxo[sContigName][SPECIES]))
    FILE.close()

def WriteReadsByFilteredContigs(dDict,dTaxo,dSeq2Quantity):
    dTemp={}
    dSize={}
    iSumReads=0
    for sContigName in dDict:
        sContigShortName=sContigName #"_".join(sContigName.split("_")[:-1])
        try:
            oCrash=dSeq2Quantity[sContigShortName]
        except KeyError:
            continue
        iReadsQuantity=dSeq2Quantity[sContigShortName]
        dSize[sContigShortName]=int(dDict[sContigName][COL_SEQSIZE])
        try:
            if sContigShortName not in dTemp[iReadsQuantity]:
                dTemp[iReadsQuantity].append(sContigShortName)
                iSumReads+=iReadsQuantity
        except KeyError:
            dTemp[iReadsQuantity]=[sContigShortName]
            iSumReads+=iReadsQuantity
    FILE=open(sOutput+"/FilteredContigs2ReadsQuantity.tsv","w")
    FILE.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("ContigId","Contig size","Reads quantity","Reads %(/Plate)","Cumulated %","Family","Genus","Species"))
    fCumulate=0.0
    for iValue in sorted(dTemp):
        for sContigName in sorted(dTemp[iValue]):
            fCumulate+=round(float(iValue)/iSumReads*100,4)
            FILE.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(sContigName,dSize[sContigName],iValue,round(float(iValue)/iSumReads*100,4),fCumulate,dTaxo[sContigName][FAMILY],dTaxo[sContigName][GENUS],dTaxo[sContigName][SPECIES]))
    FILE.close()
    
    print("Data contains {} Reads in {} Contigs".format(iSumReads,len(dDict)))
    FILE=open(sOutput+"/"+LOG,"a")
    FILE.write("Data contains {} Reads in {} Contigs\n".format(iSumReads,len(dDict)))
    FILE.close()


########################################################################
#MAIN
if __name__ == "__main__":
    print("Create Folder...")
    CreateFolder(sOutput)
    StoreLog()
    print("Load ICTV...")
    dICTV=LoadICTV(sTaxo)
    print("Load Input...")
    dBasicData=LoadData(sInput)
    print("Extract Taxo...")
    dContig2Taxo=GetTaxo(dBasicData)
    print("Load ReverseAssembly...")
    dSeq2Quantity=LoadReverseAssembly(dBasicData,sReverseAssembly)   
    print("Write Reads by Contigs...")
    WriteReadsByContigs(dBasicData,dContig2Taxo)
    print("Write Reads by Filtered Contigs...")
    WriteReadsByFilteredContigs(dBasicData,dContig2Taxo,dSeq2Quantity)
    print("Get Contigs attribut...")
    dContig2Attribut=GetContigAttribut(dBasicData)
    print("Write Stat Plates...")
    ProcessStatAllSample(sOutput,dSeq2Quantity,dContig2Taxo)
    print("Write Stat, grouped by Virus...")
    ProcessStatGroupedByVirus(sOutput,dSeq2Quantity,dContig2Taxo,dContig2Attribut)
    
	
########################################################################    
iTime2=time.time()
iDeltaTime=iTime2-iTime1
print("Script done: "+str(iDeltaTime))

