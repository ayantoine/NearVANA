# coding: utf-8
"""Python3.6"""
# compatibility: python2.7, python2.6

import time
from optparse import OptionParser

sCurrentVersionScript="v1"
iTime1=time.time()
########################################################################
'''
V1-2019/11/05
Convert fasta file to fastq file

python ConvertFastq2Fasta.py -f FASTQ
FASTQ: fastq path
'''
########################################################################
#CONSTANT

########################################################################
#Options
parser = OptionParser()
parser.add_option("-f","--fastq", dest="fastq")
parser.add_option("-o","--output", dest="output")

(options, args) = parser.parse_args()

sInput=options.fastq
if not sInput:
	sys.exit("Error : no fastq -f defined, process broken")

sOutput=options.output
if not sOutput:
	sOutput=".".join(sInput.split(".")[:-1])+".fa"
	print("Warning : no output -o defined, default : "+str(sOutput))

########################################################################
#Function 	
def ConvertFastq2Fasta(sInput,sOutput):
	FILE=open(sOutput,"w")
	iLineCounter=0
	for sLine in open(sInput):
		iLineCounter+=1
		if iLineCounter%4==1:
			sName=sLine[1:-1]
		elif iLineCounter%4==2:
			FILE.write(">{}\n{}".format(sName,sLine))
	FILE.close()

########################################################################
#MAIN
if __name__ == "__main__":
	ConvertFastq2Fasta(sInput,sOutput)
	
########################################################################    
iTime2=time.time()
iDeltaTime=iTime2-iTime1
print("Script done: "+str(iDeltaTime))

