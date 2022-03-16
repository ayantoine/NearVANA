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

python SplitSubstractedBySample.py -i INPUT

INPUT: Fastq files
'''
########################################################################
#CONSTANT
UNDERSCORE="_"

########################################################################
#Options
parser = OptionParser(conflict_handler="resolve")
parser.add_option("-i","--input", dest="input")

(options, args) = parser.parse_args()

sInput=options.input
if not sInput:
	exit("Error : no input -i defined, process broken")

########################################################################
#Function
def CreateFolder(sPath):
    try:
        os.mkdir(sPath)
    except OSError or FileExistsError:
        #print("Warning: Can not create output folder {}".format(sPath))
        pass

def ParseInput(sFastq):
    sPreviousSample=""
    sContent=""
    iIndex=0
    for sNewLine in open(sFastq):
        iIndex+=1
        sContent+=sNewLine
        if iIndex%4==1:
            sLine=sNewLine.strip()
            sSample=sLine.split(UNDERSCORE)[-1]
        if iIndex%4==0:
            if sPreviousSample!=sSample:
                try:
                    FILE.close()
                except NameError:
                    pass
                CreateFolder(sSample)
                FILE=open(sSample+"/"+sInput,"a")
            FILE.write(sContent)    
            sContent=""
    FILE.close()

########################################################################
#MAIN
if __name__ == "__main__":
    ParseInput(sInput)

########################################################################    
iTime2=time.time()
iDeltaTime=iTime2-iTime1
print("Script done: "+str(iDeltaTime))

