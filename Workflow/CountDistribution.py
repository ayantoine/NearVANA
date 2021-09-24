# coding: utf-8
"""Python3.6"""
# compatibility: python2.7, python2.6

import time
from optparse import OptionParser

sCurrentVersionScript="v1"
iTime1=time.time()
########################################################################
'''
V1-2020/01/31
Count sample occurency in file

python CountDistribution.py -i INPUT -o OUTPUT
INPUT: Demultiplexing output file
OUTPUT: Final file
'''
########################################################################
#CONSTANT
DEFAULT="UnassignedReads"

########################################################################
#Options
parser = OptionParser()
parser.add_option("-i","--input", dest="input")
parser.add_option("-o","--output", dest="output")

(options, args) = parser.parse_args()

sInput=options.input
if not sInput:
	exit("Error : no input -i defined, process broken")

sOutput=options.output
if not sOutput:
	exit("Error : no output -o defined, process broken")
	

########################################################################
#Function 	
def CountQuantity(sPath):
	dDict={}
	for sNewLine in open(sPath):
		tLine=sNewLine.split("\t")
		if len(tLine)>1:
			sSample=tLine[1]
			try:
				dDict[sSample]+=1
			except KeyError:
				dDict[sSample]=1
	return dDict

def WriteDistribution(dDict,sPath):
	FILE=open(sPath,'w')
	try:
		FILE.write(DEFAULT+"\t"+str(dDict[DEFAULT])+"\n")
	except KeyError:
		FILE.write(DEFAULT+"\t0\n")
	for sKey in sorted(dDict):
		if sKey==DEFAULT:
			continue
		try:
			FILE.write(sKey+"\t"+str(dDict[sKey])+"\n")
		except KeyError:
			FILE.write(sKey+"\t0\n")
	FILE.close()

########################################################################
#MAIN
if __name__ == "__main__":
	dSample2Quantity=CountQuantity(sInput)
	WriteDistribution(dSample2Quantity,sOutput)
	
########################################################################    
iTime2=time.time()
iDeltaTime=iTime2-iTime1
print("Script done: "+str(iDeltaTime))

