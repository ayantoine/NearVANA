# coding: utf-8
"""Python3.6"""
# compatibility: python2.7, python2.6

import time
from optparse import OptionParser

sCurrentVersionScript="v1"
iTime1=time.time()
########################################################################
'''
V2-2020/10/22
Memory consummation explode for RNAseq. Try a ligther way to do.

V1-2019/10/23
Scan two Pair-end Fastq, reorder and extract singlet.
Replace space in name by underscore "_"

python RetrievePair.py -i R1FILE -p R2FILE
R1FILE: R1.fastq
R2FILE: R2.fastq
'''
########################################################################
#CONSTANT
DEINTERLACING=".deinterlaced"

R1="R1"
R0="R0"

SPACE=" "
NOSPACE="_"

########################################################################
#Options
parser = OptionParser()
parser.add_option("-i","--input", dest="input")
parser.add_option("-p","--pair", dest="pair")

(options, args) = parser.parse_args()

sR1Fastq=options.input
if not sR1Fastq:
	exit("Error : no input -i defined, process broken")

sR2Fastq=options.pair
if not sR2Fastq:
	exit("Error : no pair -p defined, process broken")

########################################################################
#Function 	
# def GetRootName(sPath):
	# setName=set()
	# iLineCount=0
	# for sNewLine in open(sPath):
		# iLineCount+=1
		# if iLineCount%4==1:
			# sRoot=sNewLine.split(" ")[0]
			# # print(sRoot)
			# setName.add(sRoot)
	# return setName

def GetRootName(sPath):
	dDict={}
	iLineCount=0
	for sNewLine in open(sPath):
		iLineCount+=1
		if iLineCount%4==1:
			sRoot=sNewLine.split(" ")[0]
			dDict[sRoot]=True
	return dDict
	
def IntersectRootName(dBase,sPath):
	dDict={}
	iLineCount=0
	for sNewLine in open(sPath):
		iLineCount+=1
		if iLineCount%4==1:
			sRoot=sNewLine.split(" ")[0]
			try:
				oCrash=dBase[sRoot]
				dDict[sRoot]=True
				del dBase[sRoot]
			except KeyError:
				continue
				
	return dDict
	
	
def WriteFile(tListFastq,dCommon): #,setCommon):
	sFileR1=tListFastq[0]+DEINTERLACING
	sFileR2=tListFastq[1]+DEINTERLACING
	sFileR0=sFileR1.replace(R1,R0)
	
	FILER1=open(sFileR1,"w")
	FILER2=open(sFileR2,"w")
	FILER0=open(sFileR0,"w")
	
	tName=[sFileR1,sFileR2]
	tFile=[FILER1,FILER2]
	
	iSingletCount=0
	for iIndex in range(len(tListFastq)):
		sFastqFile=tListFastq[iIndex]
		oTarget=tFile[iIndex]
		iLineCounter=0
		iFileCount=0
		sSeqName=""
		sContent=""
		sInterline=""
		sQuality=""
		for sNewLine in open(sFastqFile):
			iLineCounter+=1
			if iLineCounter%4==1:
				if sSeqName!="":
					sRootName=sSeqName.split(" ")[0]
					sSeqNewName=sSeqName.replace(SPACE,NOSPACE)
					# if sRootName in setCommon:
					try:
						oCrash=dCommon[sRootName]
						iFileCount+=1
						oTarget.write(sSeqNewName+sContent+sInterline+sQuality)
					# else:
					except KeyError:
						iSingletCount+=1
						FILER0.write(sSeqNewName+sContent+sInterline+sQuality)
				sSeqName=sNewLine
				sContent=""
				sInterline=""
				sQuality=""
			elif iLineCounter%4==2:
				sContent+=sNewLine
			elif iLineCounter%4==3:
				sInterline+=sNewLine
			else:
				sQuality+=sNewLine
		if sSeqName!="":
			sRootName=sSeqName.split(" ")[0]
			sSeqNewName=sSeqName.replace(SPACE,NOSPACE)
			# if sRootName in setCommon:
			try:
				oCrash=dCommon[sRootName]
				iFileCount+=1
				oTarget.write(sSeqNewName+sContent+sInterline+sQuality)
			# else:
			except KeyError:
				iSingletCount+=1
				FILER0.write(sSeqNewName+sContent+sInterline+sQuality) 
		print(tName[iIndex]+" contains "+str(iFileCount)+" sequences")
	print(sFileR0+" contains "+str(iSingletCount)+" sequences")
			

########################################################################
#MAIN
if __name__ == "__main__":
	# setR1RootName=GetRootName(sR1Fastq)
	# print(len(setR1RootName))
	# setR2RootName=GetRootName(sR2Fastq)
	# print(len(setR2RootName))
	# setRXRootName_intersection=setR1RootName & setR2RootName
	# print(len(setRXRootName_intersection))
	dR1RootName=GetRootName(sR1Fastq)
	dRXRootName_intersection=IntersectRootName(dR1RootName,sR2Fastq)
	# WriteFile([sR1Fastq,sR2Fastq],setRXRootName_intersection)
	WriteFile([sR1Fastq,sR2Fastq],dRXRootName_intersection)
		
########################################################################    
iTime2=time.time()
iDeltaTime=iTime2-iTime1
print("Script done: "+str(iDeltaTime))

