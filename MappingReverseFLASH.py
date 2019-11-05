# coding: utf-8
"""Python3.6"""
# compatibility: python2.7, python2.6

import time
from optparse import OptionParser

sCurrentVersionScript="v1"
iTime1=time.time()
########################################################################
'''
V1-2019/11/04
Extract Sequences from Mapping and dispatch them among file.

python MappingReverseSPAdes.py -i INPUT -p PID
INPUT: path of the mapping file
PID: Project id
'''
########################################################################
#CONSTANT
CONTIG_BASENAME="Contigs_FLASH"

FLASH_ASSEMBLY_OUTPUT="FLASH.extendedFrags.fastq"

R1_TAG="1:N:0"
R2_TAG="2:N:0"

########################################################################
#Options
parser = OptionParser()
parser.add_option("-i","--input", dest="input")
parser.add_option("-p","--pid", dest="pid")

(options, args) = parser.parse_args()

sInput=options.input
if not sInput:
	exit("Error : no input -i defined, process broken")

sPID=options.pid
if not sPID:
	exit("Error : no pid -p defined, process broken")

#Half-Constant
OUTPUT_FILENAME=sPID+"_All.FLASH_reverseAssembly.tsv"
OUTPUT_FLASH_CONTIG=sPID+"_All.FLASH_contigs.fa"

########################################################################
#Function 	
def ParseAndWrite(sPath):
	FILE1=open(OUTPUT_FILENAME,"w")
	FILE2=open(OUTPUT_FLASH_CONTIG,"w")
	
	iLineCounter=0
	iContig=0
	for sNewLine in open(sPath):
		iLineCounter+=1
		if iLineCounter%4==1:
			iContig+=1
			sName=sNewLine[1:-1]
			sSample=sName.split("_")[-1]
			sNewName=CONTIG_BASENAME+"_"+str(iContig)+"_(2)"
			
			
			FILE1.write(sName+"\t"+sNewName+"\t"+sName+"\t"+sSample)
			FILE1.write(sName.replace(R1_TAG,R2_TAG)+"\t"+sNewName+"\t"+sName+"\t"+sSample)
			FILE2.write(">"+sNewName+"\n")
		elif iLineCounter%4==2:
			FILE2.write(sNewLine)
	
	FILE1.close()
	FILE2.close()


########################################################################
#MAIN
if __name__ == "__main__":
	ParseAndWrite(sInput)

########################################################################    
iTime2=time.time()
iDeltaTime=iTime2-iTime1
print("Script done: "+str(iDeltaTime))

