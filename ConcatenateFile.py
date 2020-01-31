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
Concatenate many file into one

python ConcatenateFile.py -l LISTFILE -o OUTPUT
LISTFILE: Liste of file, separated by ";"
OUTPUT: Final file
'''
########################################################################
#CONSTANT

########################################################################
#Options
parser = OptionParser()
parser.add_option("-l","--listfile", dest="listfile")
parser.add_option("-o","--output", dest="output")

(options, args) = parser.parse_args()

sListFile=options.listfile
if not sListFile:
	exit("Error : no listfile -l defined, process broken")
tListFile=sListFile.split(",")
print(tListFile)

sOutput=options.output
if not sOutput:
	exit("Error : no output -o defined, process broken")
	

########################################################################
#Function 	
def ConcatenateFile(tList,sFile):
	with open(sFile,'w') as oOutFile:
		for sPath in tList:
			with open(sPath) as oInFile:
				for sLine in oInFile:
					oOutFile.write(sLine)

########################################################################
#MAIN
if __name__ == "__main__":
	ConcatenateFile(tListFile,sOutput)
	
########################################################################    
iTime2=time.time()
iDeltaTime=iTime2-iTime1
print("Script done: "+str(iDeltaTime))

