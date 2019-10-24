# coding: utf-8
"""Python3.6"""
# compatibility: python2.7, python2.6

import time
from optparse import OptionParser
import codecs

sCurrentVersionScript="v1"
iTime1=time.time()
########################################################################
'''
V1-2019/10/24
Extract Sequences from Mapping and store them into R1 R2 or R3 data.

python MappingExtraction.py -i INPUT -p PID
INPUT: path of the mapping file
PID: Project id
'''
########################################################################
#CONSTANT
REJECT_TAG="@"
NO_MAPPING="*"
COL_SEQID=0 #! count col from 0
COL_MAPPING=2 
COL_SEQUENCE=9
COL_QUALITY=10

R1_TAG="_1:N"
R2_TAG="_2:N"
########################################################################
#Options
parser = OptionParser()
parser.add_option("-i","--input", dest="input")
parser.add_option("-p","--pid", dest="pid")

sInput=options.input
if not sInput:
	exit("Error : no input -i defined, process broken")

sPID=options.pid
if not sPID:
	exit("Error : no pid -p defined, process broken")

#Constant
R1FILE_NAME=sPID+"_R1.Substracted.fastq"
R1FILE_NAME=sPID+"_R2.Substracted.fastq"
R0FILE_NAME=sPID+"_R0.Substracted.fastq"
R3FILE_NAME=sPID+"_R1.PhiX.fastq"
R4FILE_NAME=sPID+"_R2.PhiX.fastq"

########################################################################
#Function 	
def ParseAndWrite(sPath):
	bFirst=False
	dbData=None
	
	R1FILE=open(R1FILE_NAME,"w")
	R2FILE=open(R2FILE_NAME,"w")
	R0FILE=open(R0FILE_NAME,"w")
	R3FILE=open(R3FILE_NAME,"w")
	R4FILE=open(R4FILE_NAME,"w")
	
	for sNewLine in open(sPath):
		if sNewLine[0]==REJECT_TAG:
			continue
		sLine=sNewLine.strip()
		tLine=sLine.split("\t")
		sMap=tLine[COL_MAPPING]
		sId=tLine[COL_SEQID]
		sSeq=tLine[COL_SEQUENCE]
		sQual=tLine[COL_QUALITY]
		if sMap!=NO_MAPPING:
			if R1_TAG in sId:
				R3FILE.write("@"+sId+"\n"+sSeq+"\n+\n"+sQual+"\n")
			elif R2_TAG in sId:
				R4FILE.write("@"+sId+"\n"+sSeq+"\n+\n"+sQual+"\n")
			else:
				exit("Error 80 : FATAL\nUnknown pair type: "+sSeq)
		else:
			if not bFirst:
				bFirst=True
				dbData=(sId,sSeq,sQual)
			else:
				if dbData[0].split("_")[0]==sId.split("_")[0]:
					#Same pair
					bFirst=False
					dbData=None
					if R1_TAG in sId:
						R1FILE.write("@"+sId+"\n"+sSeq+"\n+\n"+sQual+"\n")
						R2FILE.write("@"+dbData[0]+"\n"+dbData[1]+"\n+\n"+dbData[2]+"\n")
					elif R2_TAG in sId:
						R2FILE.write("@"+sId+"\n"+sSeq+"\n+\n"+sQual+"\n")
						R1FILE.write("@"+dbData[0]+"\n"+dbData[1]+"\n+\n"+dbData[2]+"\n")
					else:
						exit("Error 96 : FATAL\nUnknown pair type: "+sSeq)
				else:
					#Stored dbData is solo
					R0FILE.write("@"+dbData[0]+"\n"+dbData[1]+"\n+\n"+dbData[2]+"\n")
					bFirst=True
					dbData=(sId,sSeq,sQual)
		
	R1FILE.close()
	R2FILE.close()
	R0FILE.close()
	R3FILE.close()
	R4FILE.close()
		
########################################################################
#MAIN
if __name__ == "__main__":
	ParseAndWrite(sInput)

########################################################################    
iTime2=time.time()
iDeltaTime=iTime2-iTime1
print("Script done: "+str(iDeltaTime))

