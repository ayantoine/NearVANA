# coding: utf-8
"""Python3.6"""
# compatibility: python2.7, python2.6

import time
from optparse import OptionParser

sCurrentVersionScript="v1"
iTime1=time.time()
########################################################################
'''
V2-2019/11/21
Split a fasta by chunck of 1000 into a specified folder

python SplitFasta.py -f FOLDER -i INPUT
INPUT: Fasta file
FOLDER: Targeted folder
'''
########################################################################
#CONSTANT
CHUNCK=1000

########################################################################
#Options
parser = OptionParser()
parser.add_option("-f","--folder", dest="folder")
parser.add_option("-i","--input", dest="input")

(options, args) = parser.parse_args()

sFolder=options.folder
if not sFolder:
	exit("Error : no folder -f defined, process broken")

sInput=options.input
if not sInput:
	exit("Error : no input -i defined, process broken")

########################################################################
#Function 	
def SplitFasta(sFolder,sInput):
	iFileNumber=1
	iSeqCount=0
	FILE=open(sFolder+"/"+sInput+"."+str(iFileNumber),"w")
	sSeqName=""
	sSeqContent=""
	for sNewLine in open(sInput):
		if ">"==sNewLine[0]:
			if sSeqName!="":
				FILE.write(sSeqName+sSeqContent)
				iSeqCount+=1
				if iSeqCount==CHUNCK:
					FILE.close()
					iFileNumber+=1
					FILE=open(sFolder+"/"+sInput+"."+str(iFileNumber),"w")
				sSeqName=""
				sSeqContent=""
			sSeqName=sNewLine
		else:
			sSeqContent=sNewLine
	if sSeqName!="":
		FILE.write(sSeqName+sSeqContent)
	FILE.close()

########################################################################
#MAIN
if __name__ == "__main__":
	SplitFasta(sFolder,sInput)
	
########################################################################    
iTime2=time.time()
iDeltaTime=iTime2-iTime1
print("Script done: "+str(iDeltaTime))

