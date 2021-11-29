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
Split fasta between .keeped and .rejected based on blast result

python SplitFastaOnBlastResults.py -f FASTA -b BLAST
FASTA: fastq path
BLAST: Blast output, format 6 (tsv)
'''
########################################################################
#CONSTANT
KEEPED=".keeped"
REJECTED=".rejected"

########################################################################
#Options
parser = OptionParser()
parser.add_option("-f","--fasta", dest="fasta")
parser.add_option("-b","--blast", dest="blast")
parser.add_option("-o","--output", dest="output")

(options, args) = parser.parse_args()

sInput=options.fasta
if not sInput:
	exit("Error : no fasta -f defined, process broken")

sBlast=options.blast
if not sBlast:
	exit("Error : no blast -b defined, process broken")

sOutput=options.output
if not sOutput:
	exit("Error : no output -o defined, process broken")

########################################################################
#Function 	
def GetKeeped(sBlast):
	setKeeped=set([])
	for sNewLine in open(sBlast):
		setKeeped.add(sNewLine.split("\t")[0])
	return setKeeped
	
def SplitFasta(sFasta,setKeeped,sOutput):
	bWriteIt=False
	FILE1=open(sOutput+KEEPED,"w")
	FILE2=open(sOutput+REJECTED,"w")
	
	for sNewLine in open(sFasta):
		if ">"==sNewLine[0]:
			sName=sNewLine[1:-1]
			if sName in setKeeped:
				bWriteIt=True
			else:
				bWriteIt=False
		if bWriteIt:
			FILE1.write(sNewLine)
		else:
			FILE2.write(sNewLine)
	
	FILE1.close()
	FILE2.close()
	

########################################################################
#MAIN
if __name__ == "__main__":
	setKeepedTarget=GetKeeped(sBlast)
	SplitFasta(sInput,setKeepedTarget,sOutput)
	
########################################################################    
iTime2=time.time()
iDeltaTime=iTime2-iTime1
print("Script done: "+str(iDeltaTime))

