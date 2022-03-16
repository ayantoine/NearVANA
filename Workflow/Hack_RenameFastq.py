# coding: utf-8
"""Python3.6"""

import time
from optparse import OptionParser
import os

sCurrentVersionScript="v19"
iTime1=time.time()
########################################################################
'''
V1-2022/02/22

python Hack_RenameFastq.py -i INPUT -o OUTPUT -s SAMPLEID

INPUT: Fastq files
OUTPUT: Fastq files
SAMPLEID: SAmple id
'''
########################################################################
#CONSTANT
UNDERSCORE="_"

########################################################################
#Options
parser = OptionParser(conflict_handler="resolve")
parser.add_option("-i","--input", dest="input")
parser.add_option("-o","--output", dest="output")
parser.add_option("-s","--sample", dest="sample")

(options, args) = parser.parse_args()

sInput=options.input
if not sInput:
	exit("Error : no input -i defined, process broken")

sOutput=options.output
if not sInput:
	exit("Error : no output -o defined, process broken")

sSampleId=options.sample
if not sInput:
	exit("Error : no sample -s defined, process broken")

########################################################################
#Function
def RewriteFastq(sInput,sOutput,sSampleId):
    FILE=open(sOutput,"w")
    iIndex=0
    for sNewLine in open(sFastq):
        iIndex+=1
        sContent=sNewLine
        if iIndex%4==1:
            tContent=sContent.split(UNDERSCORE)
            tContent.append(tContent[-1])
            tContent[-2]=sSampleId
            sContent=UNDERSCORE.join(tContent)
        FILE.write(sContent)    
    FILE.close()

########################################################################
#MAIN
if __name__ == "__main__":
    RewriteFastq(sInput,sOutput,sSampleId)

########################################################################    
iTime2=time.time()
iDeltaTime=iTime2-iTime1
print("Script done: "+str(iDeltaTime))

