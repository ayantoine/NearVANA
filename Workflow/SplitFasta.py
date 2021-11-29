# coding: utf-8
"""Python3.6"""
# compatibility: python2.7, python2.6

import time
from optparse import OptionParser

sCurrentVersionScript="v3"
iTime1=time.time()
########################################################################
'''
V3-2020/01/21
Allow to define the value of CHUNCK

V2-2019/11/21
Split a fasta by chunck of 1000 into a specified folder

python SplitFasta.py -f FOLDER -i INPUT -c CHUNCK
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
parser.add_option("-c","--chunck", dest="chunck")

(options, args) = parser.parse_args()

sFolder=options.folder
if not sFolder:
	exit("Error : no folder -f defined, process broken")

sInput=options.input
if not sInput:
	exit("Error : no input -i defined, process broken")
	
sChunck=options.chunck
if not sChunck:
	iChunck=CHUNCK
	print("Warning : no chunck -c defined, default value: "+str(iChunck))
else:
	try:
		iChunck=int(sChunck)
	except ValueError:
		iChunck=CHUNCK
		print("Warning : chunck -c is not an integer, default value: "+str(iChunck))

########################################################################
#Function 	
def SplitFasta(sFolder,sInput,iChunck):
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
				if iSeqCount==iChunck:
					FILE.close()
					iSeqCount=0
					iFileNumber+=1
					FILE=open(sFolder+"/"+sInput+"."+str(iFileNumber),"w")
				sSeqName=""
				sSeqContent=""
			sSeqName=sNewLine
		else:
			sSeqContent+=sNewLine
	if sSeqName!="":
		FILE.write(sSeqName+sSeqContent)
	FILE.close()

########################################################################
#MAIN
if __name__ == "__main__":
	SplitFasta(sFolder,sInput,iChunck)
	
########################################################################    
iTime2=time.time()
iDeltaTime=iTime2-iTime1
print("Script done: "+str(iDeltaTime))

