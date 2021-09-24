# coding: utf-8
"""Python3.6"""

import time
from optparse import OptionParser
import os
import re

sCurrentVersionScript="v23"
iTime1=time.time()
########################################################################
'''
V23-2021/09/24
Adapt to new table results that fraction contigs reads among sample

V22-2021/06/18
Branching to main workflow. Master.o no more available, switch on DATA
V21-2021/04/23
Eukaryota must be subdivided into Fungi, Metazoa and Viridiplantae
V20-2021/04/20
Works for all superkingdoms
V19-2021/02/04
Retrieve Sample from Master.o to catch multiplate data.
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

python Extract4Stat.py -i INPUT -o OUTPUT -v ICTV -r REVERSE -d DATA -1 DB1 -2 DB2 -3 DB3 [-s SIZETHRESHOLD -q QUANTITYTHRESOLD -h HARDLIMIT
                        -t REJECTED -l LIMITREADS -n NOISE]
INPUT: NearVANA tsv table
OUTPUT: output folder results
ICTV: ictv taxonomy tsv table
REVERSE: ReverseAssembly file
DATA: NearVANA.data file
SIZETHRESHOLD: Minimum size to take a contig into account
QUANTITYTHRESOLD: Minimum quantity of reads to take a contig into account
HARDLIMIT: Minimum value to counter QUANTITYTHRESHOLD
REJECTED: List of taxo rejected, separated by "-". Space " " must be replace by underscore "_"
LIMITREADS: Minimum ratio of reads from a sample to assign this sample to a contig
NOISE: Part of Contigs below or equal to this value are ignored 
DB1: Genbank Family list file
DB2: Genbank Genus list file
DB3: Genbank Species list file
'''
########################################################################
#CONSTANT
BEST_HIT="Best hit"
VIRUSES="Viruses"
ALL="All"
CHAR_POINT="."

FILTER_CONTIG_SIZE=0
FILTER_READ_QUANTITY=300
FILTER_READ_REPRESENTATIVITE=1.0
FILTER_READ_HARDLIMIT=1000
FILTER_NOISE=0

INTEGER=["0","1","2","3","4","5","6","7","8","9"]
SPACE=" "
UNDERSCORE="_"

REJECTED_SEPARATOR="-"
REJECTED_TAXO=[]

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

UNCLASSIFIED="unclassified"
UNKNOWN="unknown"
UNASSIGNED="unassigned"

FAMILY="Family"
GENUS="Genus"
SPECIES="Species"

READS="Reads"
CONTIGS="Contigs"
TAXO="Taxo"
SUPERKINGDOM="Superkingdom"
SIZE="Size"
REF="Ref"
EUKARYOTA="Eukaryota"
FUNGI="Fungi"
METAZOA="Metazoa"
VIRIDIPLANTAE="Viridiplantae"
SUBCLASS=[FUNGI,METAZOA,VIRIDIPLANTAE]

EMPTY=""
SPACE=" "

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

# UNKNOWN_SIZETHRESHOLD=2000
# UNKNOWN_QUANTITYTHRESHOLD=500

FOLDER_FAMILYGLOBAL="StatByFamily"
FOLDER_GROUP="StatByGroup"
FOLDER_FAMILY="ByFamily"
FOLDER_GENUS="ByGenus"
FOLDER_SPECIES="BySpecies"
FOLDER_SAMPLE="StatBySample"
FOLDER_PLATE="StatAllSample"
FOLDER_CONTIGS="StatByContigs"
FOLDER_SAMPLE2CONTIGS="Sample2Contigs"

LOG="options.log"

DATA_COMMENT_TAG="#"
DATA_PLATE_TAG="PLATE"
DATA_PARENTHESIS_REGEX="(|)"
DATA_TABULATION="\t"

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
parser.add_option("-l","--reads_ratio_threshold", dest="reads_ratio_threshold")
parser.add_option("-n","--noise", dest="noise")
parser.add_option("-h","--hardlimit", dest="hardlimit")
parser.add_option("-1","--genbank_family", dest="genbank_family")
parser.add_option("-2","--genbank_genus", dest="genbank_genus")
parser.add_option("-3","--genbank_species", dest="genbank_species")
parser.add_option("-d","--datafile", dest="datafile")

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

sReadThreshold=options.reads_ratio_threshold
if not sReadThreshold:
    print("Warning : no reads_ratio_threshold -l defined, default value : "+str(FILTER_READ_REPRESENTATIVITE))
    sReadThreshold=str(FILTER_READ_REPRESENTATIVITE)
try:
    fReadThreshold=float(sReadThreshold)
    FILTER_READ_REPRESENTATIVITE=fReadThreshold
except ValueError:
    exit("Error : reads_ratio_threshold -l must be an integer, process broken")

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

sGBfamily=options.genbank_family
if not sGBfamily:
	exit("Error : no genbank_family -1 defined, process broken")

sGBgenus=options.genbank_genus
if not sGBgenus:
	exit("Error : no genbank_genus -2 defined, process broken")
    
sGBspecies=options.genbank_species
if not sGBspecies:
	exit("Error : no genbank_species -3 defined, process broken")

sDatafile=options.datafile
if not sDatafile:
    exit("Error : no datafile -d defined, process broken")  

########################################################################
#Function
def StoreLog():
    FILE=open(sOutput+"/"+LOG,"w")
    FILE.write("script_version: "+sCurrentVersionScript+"\n")
    FILE.write("input: "+sInput+"\n")
    FILE.write("output: "+sOutput+"\n")
    FILE.write("taxo_ictv: "+sTaxo+"\n")
    FILE.write("genbank_family: "+sGBfamily+"\n")
    FILE.write("genbank_genus: "+sGBgenus+"\n")
    FILE.write("genbank_species: "+sGBspecies+"\n")
    FILE.write("quantity_treshold: "+sQuantityThreshold+"\n")
    FILE.write("hardlimit: "+sHardLimit+"\n")
    FILE.write("size_threshold: "+sSizeThreshold+"\n")
    FILE.write("rejected_taxo: "+sRejected+"\n")
    FILE.write("reverse_assembly: "+sReverseAssembly+"\n")
    FILE.write("reads_ratio_threshold: "+sReadThreshold+"\n")
    FILE.write("reads_noise_threshold: "+sNoise+"\n")
    FILE.write("PartOf_coverThresold: "+str(COVER_THRESHOLD)+"\n")
    FILE.write("Maybe_familyThreshold: "+str(FAMILY_THRESHOLD)+"\n")
    FILE.write("Maybe_genusThreshold: "+str(GENUS_THRESHOLD)+"\n")
    FILE.write("Maybe_speciesThreshold: "+str(SPECIES_THRESHOLD)+"\n")
    FILE.close()

def CreateFolder(sPath,bAlone=False):
    if bAlone:
        tFolder=[sPath]
    else:
        tFolder=[sPath,"{}/{}".format(sPath,FOLDER_FAMILYGLOBAL),
                # "{}/{}/{}".format(sPath,FOLDER_FAMILYGLOBAL,FOLDER_GROUP),
                # "{}/{}/{}".format(sPath,FOLDER_FAMILYGLOBAL,FOLDER_SAMPLE),
                # "{}/{}".format(sPath,FOLDER_GROUP),
                # "{}/{}/{}".format(sPath,FOLDER_GROUP,FOLDER_FAMILY),
                # "{}/{}/{}".format(sPath,FOLDER_GROUP,FOLDER_GENUS),
                # "{}/{}/{}".format(sPath,FOLDER_GROUP,FOLDER_SPECIES),
                "{}/{}".format(sPath,FOLDER_SAMPLE),
                "{}/{}/{}".format(sPath,FOLDER_SAMPLE,FOLDER_FAMILY),
                "{}/{}/{}".format(sPath,FOLDER_SAMPLE,FOLDER_GENUS),
                "{}/{}/{}".format(sPath,FOLDER_SAMPLE,FOLDER_SPECIES),
                "{}/{}".format(sPath,FOLDER_PLATE),
                "{}/{}".format(sPath,FOLDER_CONTIGS),
                "{}/{}".format(sPath,FOLDER_SAMPLE2CONTIGS),
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
    for sNewLine in open(sPath):
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
    
def LoadGBdata(dDict):
    dResult={}
    dKeyword2ColId={}
    for sKey in dDict:
        for sNewLine in open(dDict[sKey]):
            sLine=sNewLine.strip()
            dResult[sLine]=sKey
    return dResult
    
def LoadData(sPath,bKingdom=True,sKingdom=""):
    dResult={}
    bHeader=True
    for sNewLine in open(sPath):
        if bHeader:
            bHeader=False
            continue
        sLine=sNewLine.strip()
        tLine=sLine.split("\t")
        bKeepIt=False
        #tLine[COL_REFSUPKINGDOM]==sKingdom
        if tLine[COL_HITRANK]==BEST_HIT and bKingdom and sKingdom+";" in tLine[COL_REFTAXONOMY] \
            and int(tLine[COL_SEQSIZE])>=FILTER_CONTIG_SIZE and int(tLine[COL_READNUMBER])>=FILTER_READ_QUANTITY:
            bKeepIt=True
        elif tLine[COL_HITRANK]==BEST_HIT and not bKingdom \
            and int(tLine[COL_SEQSIZE])>=FILTER_CONTIG_SIZE and int(tLine[COL_READNUMBER])>=FILTER_READ_QUANTITY:
            bKeepIt=True
        if bKeepIt:
            sTaxo=tLine[COL_REFTAXONOMY]
            bBreak=False
            for sRejected in REJECTED_TAXO:
                if sRejected in sTaxo:
                    bBreak=True
                    break
            if bBreak:
                continue
            sUniqSeqName=tLine[COL_QUERYID]+"_"+tLine[COL_SAMPLEID]
            dbLine=tuple(tLine)
            dResult[sUniqSeqName]=dbLine
    return dResult

def GetDateFormat(sString):
    for iIndex in range(len(sString)):
        if sString[iIndex] not in INTEGER:
            break
    if iIndex==2:
        return 2 # format: dd-mm-yyyy
    if iIndex==4:
        return 1 # format: yyyy-mm-dd

def GroupByHostDateAndLocation(dDict,sOutput):
    dSample2Date={}
    dSample2Location={}
    dSample2Host={}
    iDateFormat=None
    for iLine in dDict:
        dbData=dDict[iLine]
        sHost=dbData[COL_HOST]
        sSample=dbData[COL_SAMPLEID]
        sLocation=dbData[COL_LOCATION]
        sDate=dbData[COL_DATE]
        if iDateFormat is None:
            iDateFormat=GetDateFormat(sDate)
        if iDateFormat==1:
            sDate=sDate[:7] # format: yyyy-mm-dd
        else:
            sDate=sDate[4:] # format: dd-mm-yyyy
        dSample2Date[sSample]=sDate
        dSample2Location[sSample]=sLocation
        dSample2Host[sSample]=sHost
    dResult={}
    for sKey in dSample2Date:
        dbId=(dSample2Date[sKey],dSample2Host[sKey],dSample2Location[sKey])
        try:
            dResult[dbId].append(sKey)
        except KeyError:
            dResult[dbId]=[sKey]
    FILE=open("{}/Group2SampleId.tsv".format(sOutput),"w")
    for dbId in sorted(dResult):
        FILE.write("Grp_{}\t{}\n".format("_".join(dbId).replace("/","."),", ".join(sorted(dResult[dbId]))))
    FILE.close()
    return dResult

def ReverseTableDict(dDict):    
    dResult={}
    for sKey in dDict:
        for sValue in dDict[sKey]:
            dResult[sValue]=sKey
    return dResult

def FolderNameFormat(dbString):
    return "_".join(dbString).replace("/",".").replace(":","").replace(" ","_")

# def LoadReverseAssembly(dData,sFile):
    # dResult={}
    # for sKey in dData:
        # sContigName="_".join(sKey.split("_")[:-1])
        # sSample=sKey.split("_")[-1]        
        # try:
            # oCrash=dResult[sContigName]
        # except KeyError:
            # dResult[sContigName]={}
        # dResult[sContigName][sSample]=0
   
    # for sNewLine in open(sFile):
        # tLine=sNewLine.split("\t")
        # sContigName=tLine[1]
        # try:
            # oCrash=dResult[sContigName]
        # except KeyError:
            # continue
        # sRead=tLine[0]
        # sSample=sRead.split("_")[-1]
        # dResult[sContigName][sSample]+=1

    # if FILTER_NOISE!=0:
        # tRemoveContig=[]
        # for sContigName in dResult:
            # tRemoveSample=[]
            # #Remove reads quantity below T- sample
            # for sSample in dResult[sContigName]:
                # if dResult[sContigName][sSample]<=FILTER_NOISE:
                    # tRemoveSample.append(sSample)
            # for sSample in tRemoveSample:
                # del dResult[sContigName][sSample]
            # if len(dResult[sContigName])==0:
                # tRemoveContig.append(sContigName)
        # for sContigName in tRemoveContig:
            # del dResult[sContigName]
            
    # if FILTER_READ_REPRESENTATIVITE!=0:
        # tRemoveContig=[]
        # for sContigName in dResult:
            # #Get sum of reads
            # iSum=0
            # for sSample in dResult[sContigName]:
                # iSum+=dResult[sContigName][sSample]
                
            # tRemoveSample=[]
            # #Remove sample that have less reads that minimal value
            # for sSample in dResult[sContigName]:
                # fValue=float(dResult[sContigName][sSample])/iSum*100
                # if fValue<FILTER_READ_REPRESENTATIVITE and dResult[sContigName][sSample]<FILTER_READ_HARDLIMIT:
                    # tRemoveSample.append(sSample)
            # for sSample in tRemoveSample:
                # del dResult[sContigName][sSample]
            # if len(dResult[sContigName])==0:
                # tRemoveContig.append(sContigName)
        # for sContigName in tRemoveContig:
            # del dResult[sContigName]
            
    # return dResult

def GetTaxo(dData,dICTV,dGBtaxo):
    dResult={}
    #dTaxo={VIRUSES=dICTV,ALL=dGBtaxo}
    for sItem in dData:        
        dbLine=dData[sItem]
        sContigName=dbLine[COL_QUERYID]
        sTaxonomy=dbLine[COL_REFTAXONOMY]
        fIdentity=float(dbLine[COL_IDENTITY])
        fCover=float(dbLine[COL_QUERYCOVER])
        
        sSupKingdom=dbLine[COL_REFSUPKINGDOM]
        if sSupKingdom==VIRUSES:
            dRef=dICTV
        else:
            dRef=dGBtaxo
        
        sTaxonomy=sTaxonomy.replace("; ",";")
        tTaxonomy=sTaxonomy.split(";")
        while len(tTaxonomy[-1])==0:
            tTaxonomy=tTaxonomy[:-1]
        sFamily=None
        sGenus=None
        sSpecies=None
        for sKey in tTaxonomy:
            try:
                sTaxo=dRef[sKey]
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
        
        dResult[sContigName]={FAMILY:sFamily,GENUS:sGenus,SPECIES:sSpecies,TAXO:sTaxonomy}
    
    return dResult

def ProcessBilan(sOutput,dSeq2Sample2Quantity,tSampleList):
    iSumReads=0
    iSumContigs=0
    for sSeq in dSeq2Sample2Quantity:
        iSumContigs+=1
        for sSample in dSeq2Sample2Quantity[sSeq]:
            iSumReads+=dSeq2Sample2Quantity[sSeq][sSample]
    print("Data contains {} Reads in {} Contigs".format(iSumReads,iSumContigs))
    FILE=open(sOutput+"/"+LOG,"a")
    FILE.write("Data contains {} Reads in {} Contigs\n".format(iSumReads,iSumContigs))
    FILE.close()
    
    dSample2Count={}
    for sSeq in dSeq2Sample2Quantity:
        for sSample in dSeq2Sample2Quantity[sSeq]:
            try:
                oCrash=dSample2Count[sSample]
            except KeyError:
                dSample2Count[sSample]={CONTIGS:0,READS:0}
            dSample2Count[sSample][CONTIGS]+=1
            dSample2Count[sSample][READS]+=dSeq2Sample2Quantity[sSeq][sSample]
    
    FILE_SAMPLE=open(sOutput+"/BilanBySample.tsv","w")
    FILE_SAMPLE.write("Sample\tReads\tContigs\t%Reads(/All)\t%Contigs(/All)\n")
    for sSample in sorted(tSampleList): #sorted(dSample2Count):
        try:
            fReads=round(float(dSample2Count[sSample][READS])/iSumReads*100,2)
            fContigs=round(float(dSample2Count[sSample][CONTIGS])/iSumContigs*100,2)
            FILE_SAMPLE.write("{}\t{}\t{}\t{}\t{}\n".format(sSample,dSample2Count[sSample][READS],dSample2Count[sSample][CONTIGS],fReads,fContigs))
        except KeyError:
            FILE_SAMPLE.write("{}\t{}\t{}\t{}\t{}\n".format(sSample,0,0,0.0,0.0))
    FILE_SAMPLE.close()

def Old_ProcessBilan(sOutput,dSeq2Sample2Quantity,dGroup2Sample,tSampleList):
    iSumReads=0
    iSumContigs=0
    for sSeq in dSeq2Sample2Quantity:
        iSumContigs+=1
        for sSample in dSeq2Sample2Quantity[sSeq]:
            iSumReads+=dSeq2Sample2Quantity[sSeq][sSample]
    print("Data contains {} Reads in {} Contigs".format(iSumReads,iSumContigs))
    FILE=open(sOutput+"/"+LOG,"a")
    FILE.write("Data contains {} Reads in {} Contigs\n".format(iSumReads,iSumContigs))
    FILE.close()
    
    dSample2Group=ReverseTableDict(dGroup2Sample)
    dGroup2Count={}
    dSample2Count={}
    for sSeq in dSeq2Sample2Quantity:
        dGroup2Reads={}
        for sSample in dSeq2Sample2Quantity[sSeq]:
            dbGroup=dSample2Group[sSample]
            try:
                oCrash=dGroup2Reads[dbGroup]
            except KeyError:
                dGroup2Reads[dbGroup]=0
            try:
                oCrash=dSample2Count[sSample]
            except KeyError:
                dSample2Count[sSample]={CONTIGS:0,READS:0}
            dSample2Count[sSample][CONTIGS]+=1
            dSample2Count[sSample][READS]+=dSeq2Sample2Quantity[sSeq][sSample]
            dGroup2Reads[dbGroup]+=dSeq2Sample2Quantity[sSeq][sSample]
        for dbGroup in dGroup2Reads:
            try:
                oCrash=dGroup2Count[dbGroup]
            except KeyError:
                dGroup2Count[dbGroup]={CONTIGS:0,READS:0}
            dGroup2Count[dbGroup][CONTIGS]+=1
            dGroup2Count[dbGroup][READS]+=dGroup2Reads[dbGroup]
            
    FILE_GROUP=open(sOutput+"/BilanByGroup.tsv","w")
    FILE_GROUP.write("Group\tReads\tContigs\t%Reads(/All)\t%Contigs(/All)\n")
    for dbGroup in sorted(dGroup2Count):
        fReads=round(float(dGroup2Count[dbGroup][READS])/iSumReads*100,2)
        fContigs=round(float(dGroup2Count[dbGroup][CONTIGS])/iSumContigs*100,2)
        FILE_GROUP.write("{}\t{}\t{}\t{}\t{}\n".format(dbGroup,dGroup2Count[dbGroup][READS],dGroup2Count[dbGroup][CONTIGS],fReads,fContigs))
    FILE_GROUP.close()
    
    FILE_SAMPLE=open(sOutput+"/BilanBySample.tsv","w")
    FILE_SAMPLE.write("Group\tReads\tContigs\t%Reads(/All)\t%Contigs(/All)\n")
    for sSample in sorted(tSampleList): #sorted(dSample2Count):
        try:
            fReads=round(float(dSample2Count[sSample][READS])/iSumReads*100,2)
            fContigs=round(float(dSample2Count[sSample][CONTIGS])/iSumContigs*100,2)
            FILE_SAMPLE.write("{}\t{}\t{}\t{}\t{}\n".format(sSample,dSample2Count[sSample][READS],dSample2Count[sSample][CONTIGS],fReads,fContigs))
        except KeyError:
            FILE_SAMPLE.write("{}\t{}\t{}\t{}\t{}\n".format(sSample,0,0,0.0,0.0))
    FILE_SAMPLE.close()

def ProcessStatAllSample(sOutput,dSeq2Sample2Quantity,dContig2Taxo,dKingdom,sKingdom):
    dTaxoFamily={}
    dTaxoGenus={}
    dTaxoSpecies={}
    dTaxo2Superkingdom={}
    iSumContigs=0
    iSumReads=0
    for sSeq in dSeq2Sample2Quantity:
        iSumContigs+=1
        sFamily=dContig2Taxo[sSeq][FAMILY]
        sGenus=dContig2Taxo[sSeq][GENUS]
        sSpecies=dContig2Taxo[sSeq][SPECIES]
        if sKingdom==ALL:
            sSuperkingdom=dKingdom[sSeq]
            dTaxo2Superkingdom[sFamily]=sSuperkingdom
            dTaxo2Superkingdom[sGenus]=sSuperkingdom
            dTaxo2Superkingdom[sSpecies]=sSuperkingdom
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
        for sSample in dSeq2Sample2Quantity[sSeq]:
            dTaxoFamily[sFamily][READS]+=dSeq2Sample2Quantity[sSeq][sSample]
            dTaxoGenus[sGenus][READS]+=dSeq2Sample2Quantity[sSeq][sSample]
            dTaxoSpecies[sSpecies][READS]+=dSeq2Sample2Quantity[sSeq][sSample]
            iSumReads+=dSeq2Sample2Quantity[sSeq][sSample]
    if sKingdom!=ALL:
        FILE_FAMILY=open(sOutput+"/"+FOLDER_PLATE+"/ByFamily.tsv","w")
        FILE_FAMILY.write("Family\tReads\tContigs\t%Reads(/All)\t%Contigs(/All)\n")
        for sFamily in sorted(dTaxoFamily):
            fReads=round(float(dTaxoFamily[sFamily][READS])/iSumReads*100,2)
            fContigs=round(float(dTaxoFamily[sFamily][CONTIGS])/iSumContigs*100,2)
            FILE_FAMILY.write("{}\t{}\t{}\t{}\t{}\n".format(sFamily,dTaxoFamily[sFamily][READS],dTaxoFamily[sFamily][CONTIGS],fReads,fContigs))
        FILE_FAMILY.close()
        FILE_GENUS=open(sOutput+"/"+FOLDER_PLATE+"/ByGenus.tsv","w")
        FILE_GENUS.write("Genus\tReads\tContigs\t%Reads(/All)\t%Contigs(/All)\n")
        for sGenus in sorted(dTaxoGenus):
            fReads=round(float(dTaxoGenus[sGenus][READS])/iSumReads*100,2)
            fContigs=round(float(dTaxoGenus[sGenus][CONTIGS])/iSumContigs*100,2)
            FILE_GENUS.write("{}\t{}\t{}\t{}\t{}\n".format(sGenus,dTaxoGenus[sGenus][READS],dTaxoGenus[sGenus][CONTIGS],fReads,fContigs))
        FILE_GENUS.close()    
        FILE_SPECIES=open(sOutput+"/"+FOLDER_PLATE+"/BySpecies.tsv","w")
        FILE_SPECIES.write("Species\tReads\tContigs\t%Reads(/All)\t%Contigs(/All)\n")
        for sSpecies in sorted(dTaxoSpecies):
            fReads=round(float(dTaxoSpecies[sSpecies][READS])/iSumReads*100,2)
            fContigs=round(float(dTaxoSpecies[sSpecies][CONTIGS])/iSumContigs*100,2)
            FILE_SPECIES.write("{}\t{}\t{}\t{}\t{}\n".format(sSpecies,dTaxoSpecies[sSpecies][READS],dTaxoSpecies[sSpecies][CONTIGS],fReads,fContigs))
        FILE_SPECIES.close()
    else:
        FILE_FAMILY=open(sOutput+"/"+FOLDER_PLATE+"/ByFamily.tsv","w")
        FILE_FAMILY.write("Family\tReads\tContigs\t%Reads(/All)\t%Contigs(/All)\tSuperkingdom\n")
        for sFamily in sorted(dTaxoFamily):
            fReads=round(float(dTaxoFamily[sFamily][READS])/iSumReads*100,2)
            fContigs=round(float(dTaxoFamily[sFamily][CONTIGS])/iSumContigs*100,2)
            FILE_FAMILY.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(sFamily,dTaxoFamily[sFamily][READS],dTaxoFamily[sFamily][CONTIGS],fReads,fContigs,dTaxo2Superkingdom[sFamily]))
        FILE_FAMILY.close()
        FILE_GENUS=open(sOutput+"/"+FOLDER_PLATE+"/ByGenus.tsv","w")
        FILE_GENUS.write("Genus\tReads\tContigs\t%Reads(/All)\t%Contigs(/All)\tSuperkingdom\n")
        for sGenus in sorted(dTaxoGenus):
            fReads=round(float(dTaxoGenus[sGenus][READS])/iSumReads*100,2)
            fContigs=round(float(dTaxoGenus[sGenus][CONTIGS])/iSumContigs*100,2)
            FILE_GENUS.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(sGenus,dTaxoGenus[sGenus][READS],dTaxoGenus[sGenus][CONTIGS],fReads,fContigs,dTaxo2Superkingdom[sFamily]))
        FILE_GENUS.close()    
        FILE_SPECIES=open(sOutput+"/"+FOLDER_PLATE+"/BySpecies.tsv","w")
        FILE_SPECIES.write("Species\tReads\tContigs\t%Reads(/All)\t%Contigs(/All)\tSuperkingdom\n")
        for sSpecies in sorted(dTaxoSpecies):
            fReads=round(float(dTaxoSpecies[sSpecies][READS])/iSumReads*100,2)
            fContigs=round(float(dTaxoSpecies[sSpecies][CONTIGS])/iSumContigs*100,2)
            FILE_SPECIES.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(sSpecies,dTaxoSpecies[sSpecies][READS],dTaxoSpecies[sSpecies][CONTIGS],fReads,fContigs,dTaxo2Superkingdom[sFamily]))
        FILE_SPECIES.close()

def ProcessStatGroupedByVirus(sOutput,dSeq2Sample2Quantity,dContig2Taxo,dContig2Attribut):
    dFamily2Sample2Reads={}
    dSample2Reads={}
    dSample2Contigs={}
    dFamily2Sample2Contigs={}
    for sSeq in dSeq2Sample2Quantity:
        dSample2Bool={}
        sFamily=dContig2Taxo[sSeq][FAMILY]
        sFamily=sFamily.replace(MAYBE,"").replace(PARTOF,"")
        try:
            oCrash=dFamily2Sample2Reads[sFamily]
        except KeyError:
            dFamily2Sample2Reads[sFamily]={}
        try:
            oCrash=dFamily2Sample2Contigs[sFamily]
        except KeyError:
            dFamily2Sample2Contigs[sFamily]={}
        for sSample in dSeq2Sample2Quantity[sSeq]:
            dSample2Bool[sSample]=True
            try:
                dFamily2Sample2Reads[sFamily][sSample]+=dSeq2Sample2Quantity[sSeq][sSample]
            except KeyError:
                dFamily2Sample2Reads[sFamily][sSample]=dSeq2Sample2Quantity[sSeq][sSample]
            try:
                dSample2Reads[sSample]+=dSeq2Sample2Quantity[sSeq][sSample]
            except KeyError:
                dSample2Reads[sSample]=dSeq2Sample2Quantity[sSeq][sSample]
        for sSample in dSample2Bool:
            try:
                dSample2Contigs[sSample]+=1
            except KeyError:
                dSample2Contigs[sSample]=1
            try:
                dFamily2Sample2Contigs[sFamily][sSample].append(sSeq)
            except KeyError:
                dFamily2Sample2Contigs[sFamily][sSample]=[sSeq]
        
    for sFamily in dFamily2Sample2Contigs:
        # print(sOutput+"/"+FOLDER_FAMILYGLOBAL+"/"+FOLDER_SAMPLE+"/"+sFamily.replace(" ","_")+".tsv")
        # FILE_SAMPLE=open(sOutput+"/"+FOLDER_FAMILYGLOBAL+"/"+FOLDER_SAMPLE+"/"+sFamily.replace(" ","_")+".tsv","w")
        FILE_SAMPLE=open(sOutput+"/"+FOLDER_FAMILYGLOBAL+"/"+sFamily.replace(" ","_")+".tsv","w")
        FILE_SAMPLE.write("Sample\tReads\tContigs\t%Reads(/Sample)\t%Contigs(/Sample)\tContigs\tContigsSize\tRefId\n")
        for sSample in sorted(dFamily2Sample2Contigs[sFamily]):
            fReads=round(float(dFamily2Sample2Reads[sFamily][sSample])/dSample2Reads[sSample]*100,2)
            fContigs=round(float(len(dFamily2Sample2Contigs[sFamily][sSample]))/dSample2Contigs[sSample]*100,2)
            FILE_SAMPLE.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                                        sSample,dFamily2Sample2Reads[sFamily][sSample],
                                        len(dFamily2Sample2Contigs[sFamily][sSample]),fReads,fContigs,
                                        ";".join(dFamily2Sample2Contigs[sFamily][sSample]),
                                        ";".join([dContig2Attribut[X][SIZE] for X in dFamily2Sample2Contigs[sFamily][sSample]]),
                                        ";".join([dContig2Attribut[X][REF] for X in dFamily2Sample2Contigs[sFamily][sSample]])
                                        ))
        FILE_SAMPLE.close()

def Old_ProcessStatGroupedByVirus(sOutput,dSeq2Sample2Quantity,dContig2Taxo,dGroup2Sample,dContig2Attribut):
    dSample2Group=ReverseTableDict(dGroup2Sample)
    dFamily2Group2Reads={}
    dFamily2Sample2Reads={}
    dGroup2Reads={}
    dSample2Reads={}
    dGroup2Contigs={}
    dSample2Contigs={}
    dFamily2Group2Contigs={}
    dFamily2Sample2Contigs={}
    for sSeq in dSeq2Sample2Quantity:
        dGroup2Bool={}
        dSample2Bool={}
        sFamily=dContig2Taxo[sSeq][FAMILY]
        sFamily=sFamily.replace(MAYBE,"").replace(PARTOF,"")
        try:
            oCrash=dFamily2Group2Reads[sFamily]
        except KeyError:
            dFamily2Group2Reads[sFamily]={}
        try:
            oCrash=dFamily2Sample2Reads[sFamily]
        except KeyError:
            dFamily2Sample2Reads[sFamily]={}
        try:
            oCrash=dFamily2Group2Contigs[sFamily]
        except KeyError:
            dFamily2Group2Contigs[sFamily]={}
        try:
            oCrash=dFamily2Sample2Contigs[sFamily]
        except KeyError:
            dFamily2Sample2Contigs[sFamily]={}
        for sSample in dSeq2Sample2Quantity[sSeq]:
            dbGroup=dSample2Group[sSample]
            dGroup2Bool[dbGroup]=True
            dSample2Bool[sSample]=True
            try:
                dFamily2Group2Reads[sFamily][dbGroup]+=dSeq2Sample2Quantity[sSeq][sSample]
            except KeyError:
                dFamily2Group2Reads[sFamily][dbGroup]=dSeq2Sample2Quantity[sSeq][sSample]
            try:
                dFamily2Sample2Reads[sFamily][sSample]+=dSeq2Sample2Quantity[sSeq][sSample]
            except KeyError:
                dFamily2Sample2Reads[sFamily][sSample]=dSeq2Sample2Quantity[sSeq][sSample]
            try:
                dGroup2Reads[dbGroup]+=dSeq2Sample2Quantity[sSeq][sSample]
            except KeyError:
                dGroup2Reads[dbGroup]=dSeq2Sample2Quantity[sSeq][sSample]
            try:
                dSample2Reads[sSample]+=dSeq2Sample2Quantity[sSeq][sSample]
            except KeyError:
                dSample2Reads[sSample]=dSeq2Sample2Quantity[sSeq][sSample]
        for dbGroup in dGroup2Bool:
            try:
                dGroup2Contigs[dbGroup]+=1
            except KeyError:
                dGroup2Contigs[dbGroup]=1
            try:
                dFamily2Group2Contigs[sFamily][dbGroup]+=1
            except KeyError:
                dFamily2Group2Contigs[sFamily][dbGroup]=1
        for sSample in dSample2Bool:
            try:
                dSample2Contigs[sSample]+=1
            except KeyError:
                dSample2Contigs[sSample]=1
            try:
                dFamily2Sample2Contigs[sFamily][sSample].append(sSeq)
            except KeyError:
                dFamily2Sample2Contigs[sFamily][sSample]=[sSeq]
    
    for sFamily in dFamily2Group2Contigs:
        FILE_GROUP=open(sOutput+"/"+FOLDER_FAMILYGLOBAL+"/"+FOLDER_GROUP+"/"+sFamily.replace(" ","_")+".tsv","w")
        FILE_GROUP.write("Groups\tReads\tContigs\t%Reads(/Group)\t%Contigs(/Group)\t#Samples\n")
        for dbGroup in sorted(dFamily2Group2Contigs[sFamily]):
            iGroup=len(dGroup2Sample[dbGroup])
            fReads=round(float(dFamily2Group2Reads[sFamily][dbGroup])/dGroup2Reads[dbGroup]*100,2)
            fContigs=round(float(dFamily2Group2Contigs[sFamily][dbGroup])/dGroup2Contigs[dbGroup]*100,2)
            FILE_GROUP.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(
                                        dbGroup,dFamily2Group2Reads[sFamily][dbGroup],
                                        dFamily2Group2Contigs[sFamily][dbGroup],fReads,fContigs,iGroup
                                        ))
        FILE_GROUP.close()
        
    for sFamily in dFamily2Sample2Contigs:
        # print(sOutput+"/"+FOLDER_FAMILYGLOBAL+"/"+FOLDER_SAMPLE+"/"+sFamily.replace(" ","_")+".tsv")
        FILE_SAMPLE=open(sOutput+"/"+FOLDER_FAMILYGLOBAL+"/"+FOLDER_SAMPLE+"/"+sFamily.replace(" ","_")+".tsv","w")
        FILE_SAMPLE.write("Sample\tReads\tContigs\t%Reads(/Sample)\t%Contigs(/Sample)\tContigs\tContigsSize\tRefId\n")
        for sSample in sorted(dFamily2Sample2Contigs[sFamily]):
            fReads=round(float(dFamily2Sample2Reads[sFamily][sSample])/dSample2Reads[sSample]*100,2)
            fContigs=round(float(len(dFamily2Sample2Contigs[sFamily][sSample]))/dSample2Contigs[sSample]*100,2)
            FILE_SAMPLE.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                                        sSample,dFamily2Sample2Reads[sFamily][sSample],
                                        len(dFamily2Sample2Contigs[sFamily][sSample]),fReads,fContigs,
                                        ";".join(dFamily2Sample2Contigs[sFamily][sSample]),
                                        ";".join([dContig2Attribut[X][SIZE] for X in dFamily2Sample2Contigs[sFamily][sSample]]),
                                        ";".join([dContig2Attribut[X][REF] for X in dFamily2Sample2Contigs[sFamily][sSample]])
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

def ProcessStatBySample(sOutput,dSeq2Sample2Quantity,dContig2Taxo,tSampleList,dKingdom,sKingdom):
    dTaxo2Sample2Data={FAMILY:{},GENUS:{},SPECIES:{}}
    dSample2SumContigs={}
    dSample2SumReads={}
    dTaxo2Superkingdom={}
    for sSeq in dSeq2Sample2Quantity:
        dSample2Bool={}
        sFamily=dContig2Taxo[sSeq][FAMILY]
        sGenus=dContig2Taxo[sSeq][GENUS]
        sSpecies=dContig2Taxo[sSeq][SPECIES]
        if sKingdom==ALL:
            dTaxo2Superkingdom[sFamily]=dKingdom[sSeq]
            dTaxo2Superkingdom[sGenus]=dKingdom[sSeq]
            dTaxo2Superkingdom[sSpecies]=dKingdom[sSeq]
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
    
    for sSample in tSampleList:
        dTarget={FAMILY:FOLDER_FAMILY,GENUS:FOLDER_GENUS,SPECIES:FOLDER_SPECIES}
        try:
            iSumContigs=dSample2SumContigs[sSample]
            iSumReads=dSample2SumReads[sSample]
        except KeyError:
            pass
        if sKingdom!=ALL:
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
        else:
            for sTaxo in dTarget:
                FILE_SAMPLE=open(sOutput+"/"+FOLDER_SAMPLE+"/"+dTarget[sTaxo]+"/"+sSample+".tsv","w")
                FILE_SAMPLE.write(sTaxo+"\tReads\tContigs\t%Reads(/Sample)\t%Contigs(/Sample)\tSuperkingdom\n")
                try:
                    for sItem in sorted(dTaxo2Sample2Data[sTaxo][sSample]):
                        iSampleContigs=dTaxo2Sample2Data[sTaxo][sSample][sItem][CONTIGS]
                        iSampleReads=dTaxo2Sample2Data[sTaxo][sSample][sItem][READS]
                        fContigs=round(float(iSampleContigs)/iSumContigs*100,2)
                        fReads=round(float(iSampleReads)/iSumReads*100,2)
                        FILE_SAMPLE.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(sItem,iSampleReads,iSampleContigs,fReads,fContigs,dTaxo2Superkingdom[sItem]))
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

def WriteReadsByContigs(dDict,dTaxo,sOutput,dKingdom,sKingdom):
    dTemp={}
    dSize={}
    iSumReads=0
    for sContigName in dDict:
        iReadsQuantity=int(sContigName.split("(")[-1].split(")")[0])
        sContigShortName="_".join(sContigName.split("_")[:-1])
        dSize[sContigShortName]=int(dDict[sContigName][COL_SEQSIZE])
        try:
            if sContigShortName not in dTemp[iReadsQuantity]:
                dTemp[iReadsQuantity].append(sContigShortName)
                iSumReads+=iReadsQuantity
        except KeyError:
            dTemp[iReadsQuantity]=[sContigShortName]
            iSumReads+=iReadsQuantity
    FILE=open(sOutput+"/Contigs2ReadsQuantity.tsv","w")
    if sKingdom!=ALL:
        FILE.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("ContigId","Contig size","Reads quantity","Reads %(/Plate)","Cumulated %","Family","Genus","Species","Taxonomy"))
    else:
        FILE.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("ContigId","Contig size","Reads quantity","Reads %(/Plate)","Cumulated %","Superkingdom","Family","Genus","Species","Taxonomy"))
    fCumulate=0.0
    for iValue in sorted(dTemp):
        for sContigName in sorted(dTemp[iValue]):
            fCumulate+=round(float(iValue)/iSumReads*100,4)
            if sKingdom!=ALL:
                FILE.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(sContigName,dSize[sContigName],iValue,round(float(iValue)/iSumReads*100,4),fCumulate,dTaxo[sContigName][FAMILY],dTaxo[sContigName][GENUS],dTaxo[sContigName][SPECIES],dTaxo[sContigName][TAXO]))
            else:
                FILE.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(sContigName,dSize[sContigName],iValue,round(float(iValue)/iSumReads*100,4),fCumulate,dKingdom[sContigName],dTaxo[sContigName][FAMILY],dTaxo[sContigName][GENUS],dTaxo[sContigName][SPECIES],dTaxo[sContigName][TAXO]))
    FILE.close()

def WriteReadsByFilteredContigs(dDict,dTaxo,dSeq2Sample2Quantity,sOutput,dKingdom,sKingdom):
    dTemp={}
    dSize={}
    iSumReads=0
    for sContigName in dDict:
        sContigShortName="_".join(sContigName.split("_")[:-1])
        #iReadsQuantity=int(sContigName.split("(")[-1].split(")")[0])
        iReadsQuantity=0
        try:
            oCrash=dSeq2Sample2Quantity[sContigShortName]
        except KeyError:
            continue
        for sSample in dSeq2Sample2Quantity[sContigShortName]:
            iReadsQuantity+=dSeq2Sample2Quantity[sContigShortName][sSample]
        dSize[sContigShortName]=int(dDict[sContigName][COL_SEQSIZE])
        try:
            if sContigShortName not in dTemp[iReadsQuantity]:
                dTemp[iReadsQuantity].append(sContigShortName)
                iSumReads+=iReadsQuantity
        except KeyError:
            dTemp[iReadsQuantity]=[sContigShortName]
            iSumReads+=iReadsQuantity
    FILE=open(sOutput+"/FilteredContigs2ReadsQuantity.tsv","w")
    if sKingdom!=ALL:
        FILE.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("ContigId","Contig size","Reads quantity","Reads %(/Plate)","Cumulated %","Family","Genus","Species","Taxonomy"))
    else:
        FILE.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("ContigId","Contig size","Reads quantity","Reads %(/Plate)","Cumulated %","Superkingdom","Family","Genus","Species","Taxonomy"))
    fCumulate=0.0
    for iValue in sorted(dTemp):
        for sContigName in sorted(dTemp[iValue]):
            fCumulate+=round(float(iValue)/iSumReads*100,4)
            if sKingdom!=ALL:
                FILE.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(sContigName,dSize[sContigName],iValue,round(float(iValue)/iSumReads*100,4),fCumulate,dTaxo[sContigName][FAMILY],dTaxo[sContigName][GENUS],dTaxo[sContigName][SPECIES],dTaxo[sContigName][TAXO]))
            else:
                FILE.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(sContigName,dSize[sContigName],iValue,round(float(iValue)/iSumReads*100,4),fCumulate,dKingdom[sContigName],dTaxo[sContigName][FAMILY],dTaxo[sContigName][GENUS],dTaxo[sContigName][SPECIES],dTaxo[sContigName][TAXO]))
    FILE.close()

def WriteFilteredContigsData(sOutput,dDict):
    for sContigName in dDict:
        FILE=open(sOutput+"/"+FOLDER_CONTIGS+"/"+sContigName+".tsv","w")
        FILE.write("{}\t{}\t{}\n".format("SampleId","Reads","Reads %"))
        iSum=0
        for sSampleId in dDict[sContigName]:
            iSum+=dDict[sContigName][sSampleId]
        for sSampleId in sorted(dDict[sContigName]):
            FILE.write("{}\t{}\t{}\n".format(sSampleId,dDict[sContigName][sSampleId],round(float(dDict[sContigName][sSampleId])/iSum*100,4)))
        FILE.close()


def ProcessSample2Contigs(sOutput,dSeq2Sample2Quantity,dContig2Taxo,tSampleList,dSize,dKingdom,sKingdom):
    #Retrieve data
    dResult={}
    dContigs2FilteredQuantity={}
    for sSeq in dSeq2Sample2Quantity:
        try:
            oCrash=dContigs2FilteredQuantity[sSeq]
        except KeyError:
            dContigs2FilteredQuantity[sSeq]=0
        for sSample in dSeq2Sample2Quantity[sSeq]:
            try:
                oCrash=dResult[sSample]
            except KeyError:
                dResult[sSample]={READS:0}
            try:
                oCrash=dResult[sSeq]
            except KeyError:
                dResult[sSample][sSeq]={READS:dSeq2Sample2Quantity[sSeq][sSample],
                                        TAXO:{FAMILY:dContig2Taxo[sSeq][FAMILY],
                                        GENUS:dContig2Taxo[sSeq][GENUS],
                                        SPECIES:dContig2Taxo[sSeq][SPECIES]},
                                        }
            dResult[sSample][READS]+=dSeq2Sample2Quantity[sSeq][sSample]
            dContigs2FilteredQuantity[sSeq]+=dSeq2Sample2Quantity[sSeq][sSample]
    #Write data
    for sSample in tSampleList: #dResult:
        FILE_SAMPLE=open(sOutput+"/"+FOLDER_SAMPLE2CONTIGS+"/"+sSample+".tsv","w")
        if sKingdom==ALL:
            FILE_SAMPLE.write("ContigId\tContigSize\tReads\t%Reads(/Sample)\t%Reads(/Contigs)\tSuperkingdom\tFamily\tGenus\tSpecies\n")
        else:
            FILE_SAMPLE.write("ContigId\tContigSize\tReads\t%Reads(/Sample)\t%Reads(/Contigs)\tFamily\tGenus\tSpecies\n")
        try:
            for sSeq in sorted(dResult[sSample]):
                if sSeq==READS:
                    continue
                iContigSize=dSize[sSeq]
                iCurrentReads=dResult[sSample][sSeq][READS]
                iSampleReads=dResult[sSample][READS]
                iContigReads=dContigs2FilteredQuantity[sSeq]
                fReadsBySample=round(float(iCurrentReads)/iSampleReads*100,2)
                fReadsByContig=round(float(iCurrentReads)/iContigReads*100,2)
                sFamily=dContig2Taxo[sSeq][FAMILY]
                sGenus=dContig2Taxo[sSeq][GENUS]
                sSpecies=dContig2Taxo[sSeq][SPECIES]
                if sKingdom==ALL:    
                    sSuperkingdom=dKingdom[sSeq]
                    FILE_SAMPLE.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(sSeq,
                                                                    iContigSize,
                                                                    iCurrentReads,
                                                                    fReadsBySample,
                                                                    fReadsByContig,
                                                                    sSuperkingdom,
                                                                    sFamily,
                                                                    sGenus,
                                                                    sSpecies))
                else:
                    FILE_SAMPLE.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(sSeq,
                                                                    iContigSize,
                                                                    iCurrentReads,
                                                                    fReadsBySample,
                                                                    fReadsByContig,
                                                                    sFamily,
                                                                    sGenus,
                                                                    sSpecies))
        except KeyError:
            pass
        FILE_SAMPLE.close()
 
def GetSampleList(sPath):
    print(sPath)
    dTag2File={}
    tTag=[]
    for sNewLine in open(sPath):
        sLine=sNewLine.strip()
        if len(sLine)==0:
            continue
        if sLine[0]==DATA_COMMENT_TAG:
            continue
        if DATA_PLATE_TAG in sLine:
            sPart2=sLine.split(EQUAL)[-1]
            sContent=EMPTY.join(re.split(DATA_PARENTHESIS_REGEX,sPart2))
            tTag=sContent.split(SPACE)
            continue
        for sTag in tTag:
            if sTag in sLine:
                sPart2=sLine.split(EQUAL)[-1]
                sContent=EMPTY.join(re.split(DATA_PARENTHESIS_REGEX,sLine))
                tContent=sContent.split(SPACE)
                sFile=tContent[-1]
                dTag2File[sTag]=sFile
    tResult=[]
    for sTag in sorted(dTag2File):
        for sNewLine in open(dTag2File[sTag]):
            sLine=sNewLine.strip()
            if len(sLine)==0:
                continue
            tLine=sLine.split(DATA_TABULATION)
            sSampleId=sTag+tLine[0]
            tResult.append(sSampleId)
    return tResult

def GetSupkingdomList(sPath,tList):
    dResult={}
    bHeader=True
    for sNewLine in open(sPath):
        if bHeader:
            bHeader=False
            continue
        sLine=sNewLine.strip()
        tLine=sLine.split("\t")
        if tLine[COL_HITRANK]==BEST_HIT \
            and int(tLine[COL_SEQSIZE])>=FILTER_CONTIG_SIZE and int(tLine[COL_READNUMBER])>=FILTER_READ_QUANTITY:
            sAddKingdom=tLine[COL_REFSUPKINGDOM]
            if sAddKingdom==EUKARYOTA:
                sTaxo=tLine[COL_REFTAXONOMY]
                sAddKingdom=CheckSubClassTaxo(sTaxo)
            dResult[sAddKingdom]=True
    for sKey in sorted(dResult):
        if sKey!=CHAR_POINT:
            tList.append(sKey)
    return tList

def GetSize(dDict):
    dResult={}
    for sItem in dDict:        
        dbLine=dDict[sItem]
        sContigName=dbLine[COL_QUERYID]
        iSize=int(dbLine[COL_SEQSIZE])
        dResult[sContigName]=iSize
    return dResult
            
def GetKingdom(dDict):
    dResult={}
    for sItem in dDict:        
        dbLine=dDict[sItem]
        sContigName=dbLine[COL_QUERYID]
        sKingdom=dbLine[COL_REFSUPKINGDOM]
        if sKingdom==EUKARYOTA:
            sTaxo=dbLine[COL_REFTAXONOMY]
            sKingdom=CheckSubClassTaxo(sTaxo)
        dResult[sContigName]=sKingdom
    return dResult

def CheckSubClassTaxo(sString):
    bBreak=False
    for sSubTarget in SUBCLASS:
        if sSubTarget in sString:
            bBreak=True
            break
    if bBreak:
        return sSubTarget
    else:
        sNewString=sString.replace("; ",";")
        tNewString=sNewString.split(";")
        return tNewString[1]

def RetrieveReverseAssembly(dData):
    dResult={}
    #Retrieve basic data
    for sKey in dData:
        sContigName=dData[sKey][COL_QUERYID]
        try:
            oCrash=dResult[sContigName]
        except KeyError:
            dResult[sContigName]={}
        sSample=dData[sContigName][COL_SAMPLEID]
        iReads=int(dData[sContigName][COL_READNUMBER])
        try:
            dResult[sContigName][sSample]+=iReads
        except KeyError:
            dResult[sContigName][sSample]=iReads
    
    #Filter data below Noise value
    if FILTER_NOISE!=0:
        tRemoveContig=[]
        for sContigName in dResult:
            tRemoveSample=[]
            #Remove reads quantity below T- sample
            for sSample in dResult[sContigName]:
                if dResult[sContigName][sSample]<=FILTER_NOISE:
                    tRemoveSample.append(sSample)
            for sSample in tRemoveSample:
                del dResult[sContigName][sSample]
            if len(dResult[sContigName])==0:
                tRemoveContig.append(sContigName)
        for sContigName in tRemoveContig:
            del dResult[sContigName]
    
    #Filter sample association if reads below minimum representativite
    if FILTER_READ_REPRESENTATIVITE!=0:
        tRemoveContig=[]
        for sContigName in dResult:
            #Get sum of reads
            iSum=0
            for sSample in dResult[sContigName]:
                iSum+=dResult[sContigName][sSample]
                
            tRemoveSample=[]
            #Remove sample that have less reads that minimal value
            for sSample in dResult[sContigName]:
                fValue=float(dResult[sContigName][sSample])/iSum*100
                if fValue<FILTER_READ_REPRESENTATIVITE and dResult[sContigName][sSample]<FILTER_READ_HARDLIMIT:
                    tRemoveSample.append(sSample)
            for sSample in tRemoveSample:
                del dResult[sContigName][sSample]
            if len(dResult[sContigName])==0:
                tRemoveContig.append(sContigName)
        for sContigName in tRemoveContig:
            del dResult[sContigName]
            
    return dResult

########################################################################
#MAIN
if __name__ == "__main__":
    print("Create Folder...")
    CreateFolder(sOutput,bAlone=True)
    StoreLog()
    print("Load ICTV...")
    dICTV=LoadICTV(sTaxo)
    print("Load Genbank Taxonomy...")
    dGBtaxo=LoadGBdata({FAMILY:sGBfamily,GENUS:sGBgenus,SPECIES:sGBspecies})
    print("Retrieve all sample Id")
    tSampleList=GetSampleList(sDatafile)   
    print("Retrieve all superkingdom")
    tSupKingdom=GetSupkingdomList(sInput,[ALL])
    print(tSupKingdom)
    for sKingdom in tSupKingdom:
        print("---------->Working on {}<----------".format(sKingdom))
        print("Create dedicated Folder...")
        sDedicatedOutput=sOutput+"/"+sKingdom.replace(SPACE,UNDERSCORE)
        CreateFolder(sDedicatedOutput)
        print("Load Input...")
        if sKingdom==ALL:
            dBasicData=LoadData(sInput,bKingdom=False) 
            print("Retrieve superkingdom...")
            dSeq2Kingdom=GetKingdom(dBasicData)   
            dSeq2Sample2Quantity=RetrieveReverseAssembly(dBasicData)
        else:
            dBasicData=LoadData(sInput,True,sKingdom) 
            dSeq2Kingdom={}
            dSeq2Sample2Quantity=RetrieveReverseAssembly(dBasicData)
        print("Extract Taxo...")
        dContig2Taxo=GetTaxo(dBasicData,dICTV,dGBtaxo)       
        # print("Load ReverseAssembly...")
        # dSeq2Sample2Quantity=LoadReverseAssembly(dBasicData,sReverseAssembly)   
        print("Retrieve contig size...")
        dSeq2Size=GetSize(dBasicData)   
        print("Write Filtered Contigs data...")
        WriteFilteredContigsData(sDedicatedOutput,dSeq2Sample2Quantity)
        print("Write Reads by Contigs...")
        WriteReadsByContigs(dBasicData,dContig2Taxo,sDedicatedOutput,dSeq2Kingdom,sKingdom)
        print("Write Reads by Filtered Contigs...")
        WriteReadsByFilteredContigs(dBasicData,dContig2Taxo,dSeq2Sample2Quantity,sDedicatedOutput,dSeq2Kingdom,sKingdom)
        print("Get Contigs attribut...")
        dContig2Attribut=GetContigAttribut(dBasicData)
        print("Write Bilan...")
        ProcessBilan(sDedicatedOutput,dSeq2Sample2Quantity,tSampleList)
        print("Write Stat All Sample merged...")
        ProcessStatAllSample(sDedicatedOutput,dSeq2Sample2Quantity,dContig2Taxo,dSeq2Kingdom,sKingdom)
        print("Write Stat, grouped by Virus...")
        ProcessStatGroupedByVirus(sDedicatedOutput,dSeq2Sample2Quantity,dContig2Taxo,dContig2Attribut)
        print("Write Stat by Sample...")
        ProcessStatBySample(sDedicatedOutput,dSeq2Sample2Quantity,dContig2Taxo,tSampleList,dSeq2Kingdom,sKingdom)
        print("Write Sample to Contigs...")
        ProcessSample2Contigs(sDedicatedOutput,dSeq2Sample2Quantity,dContig2Taxo,tSampleList,dSeq2Size,dSeq2Kingdom,sKingdom)
    
	
########################################################################    
iTime2=time.time()
iDeltaTime=iTime2-iTime1
print("Script done: "+str(iDeltaTime))

