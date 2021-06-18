# coding: utf-8
"""Python3.6"""

import time
from optparse import OptionParser
import os

sCurrentVersionScript="v1"
iTime1=time.time()
########################################################################
'''
V1-2019/09/10

python ExtractTaxonomyFromGenbankData.py.py -1 NODES -2 NAMES -o OUTPUT -t TARGET
NODES: nodes.dmp
NAMES: names.dmp
TARGET: taxonomy level to target
OUTPUT: output file
'''
########################################################################
#CONSTANT
TABULATION="\t"
SCIENTIFIC_NAME="scientific name"

DEBUG_NUMBER=10000000

########################################################################
#Options
parser = OptionParser(conflict_handler="resolve")
parser.add_option("-1","--nodesfile", dest="nodesfile")
parser.add_option("-2","--namesfile", dest="namesfile")
parser.add_option("-o","--output", dest="output")
parser.add_option("-t","--target", dest="target")

(options, args) = parser.parse_args()

sNodeFile=options.nodesfile
if not sNodeFile:
    exit("Error : no nodesfile -1 defined, process broken")

sNameFile=options.namesfile
if not sNameFile:
    exit("Error : no namesfile -2 defined, process broken")
    
sTarget=options.target
if not sTarget:
    exit("Error : no target -t defined, process broken")

sOutput=options.output
if not sOutput:
    exit("Error : no output -o defined, process broken")

########################################################################
#Function
def TellMe(iValue):
    iValue+=1
    if iValue%DEBUG_NUMBER==0:
        print("{}...".format(iValue))
    return iValue

def ExtractTaxId(sString,sFile):
    dResult={}
    iCounter=0
    for sNewLine in open(sFile):
        iCounter=TellMe(iCounter)
        if sString in sNewLine:
            dResult[sNewLine.split(TABULATION)[0]]=True
    return dResult

def ExtractTaxonomyPart(dDict,sFile):
    dResult={}
    iCounter=0
    for sNewLine in open(sFile):
        iCounter=TellMe(iCounter)
        if SCIENTIFIC_NAME in sNewLine:
            tLine=sNewLine.split(TABULATION)
            sTaxId=tLine[0]
            sTaxName=tLine[2]
            try:
                oCrash=dDict[sTaxId]
                dResult[sTaxName]=True
            except KeyError:
                continue
    return dResult

def WriteOutput(dDict,sFile):
    FILE=open(sFile,"w")
    for sKey in sorted(dDict):
        FILE.write("{}\n".format(sKey))
    FILE.close()

########################################################################
#MAIN
if __name__ == "__main__":
    dTaxId=ExtractTaxId(sTarget,sNodeFile)
    dTarget=ExtractTaxonomyPart(dTaxId,sNameFile)
    WriteOutput(dTarget,sOutput)
    
########################################################################    
iTime2=time.time()
iDeltaTime=iTime2-iTime1
print("Script done: "+str(iDeltaTime))

