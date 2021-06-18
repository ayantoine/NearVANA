# coding: utf-8
"""Python3.6"""

import time
from optparse import OptionParser
import requests
import re

sCurrentVersionScript="v1"
iTime1=time.time()
########################################################################
'''
V1-2020/02/17
Extract viral data from ViralZone

python GetMeanViralLength.py -f FOLDER
FOLDER: output folder where to store the table
'''
########################################################################
#CONSTANT
OUTPUT="ViralFamily2MinLen.tsv"

VIRALZONE_SITE="https://viralzone.expasy.org"
VIRALZONE_TABLE=["238","239","246","243","245","240"]

ENDLINE="\n"
TABULATION="\t"

INTEGER=["0","1","2","3","4","5","6","7","8","9"]

FAMILY_GENUS="Family/genus"
SPECIES="species"
END_TABLE="</TABLE>"
HREF="href"
SPACE=" "
SPACE_B=" b"
EMPTY=""
POINT="."
COMA=","
VIRIDAE="viridae"
VIRUS="virus"
CIRCULAR="Circular"
LINEAR="Linear"

CONVERT={"kb":1000,"gb":1000000}

REGEX_BALISE="<|>"

########################################################################
#Options
parser = OptionParser()
parser.add_option("-f","--folder", dest="folder")

(options, args) = parser.parse_args()

sFolder=options.folder
# if not sFolder:
    # exit("Error : no folder -f defined, process broken")

########################################################################
#Function
def HaveNumerical(sString):
    for sInteger in INTEGER:
        if sInteger in sString:
            return True
    return False

def ParseViralZoneTable():
    dResult={}
    for sSuffix in VIRALZONE_TABLE:
        # print("---------------"+sSuffix+"---------------")
        oPage=requests.get(VIRALZONE_SITE+"/"+sSuffix)
        bOk=False
        bFirst=True
        sTarget=None
        for sNewLine in oPage.text.split(ENDLINE):
            sLine=sNewLine.strip()
            if len(sLine)==0:
                continue
            if FAMILY_GENUS in sLine:
                bOk=True
                continue
            if bOk and bFirst:
                bFirst=False
                continue
            if bOk and END_TABLE in sLine:
                bOk=False
            if bOk:
                if VIRIDAE in sLine or VIRUS in sLine:
                    tTemp=re.split(REGEX_BALISE,sLine)
                    for sItem in tTemp:
                        if VIRIDAE in sItem or VIRUS in sItem:
                            sTarget=sItem
                            continue
                else:
                    if SPECIES in sLine:
                        continue
                    # print(sLine)
                    tTemp=re.split(REGEX_BALISE,sLine)
                    iOccur=0
                    for sItem in tTemp:
                        if HaveNumerical(sItem):
                            if CIRCULAR in sItem or LINEAR in sItem:
                                continue
                            iOccur+=1
                            # print(iOccur,sItem)
                            if iOccur==1:
                                oMin=sItem.replace(SPACE_B,EMPTY).replace(SPACE,EMPTY).replace(COMA,POINT).lower()
                            elif iOccur==2:
                                oMax=sItem.replace(SPACE_B,EMPTY).replace(SPACE,EMPTY).replace(COMA,POINT).lower()
                    for sKey in CONVERT:
                        if sKey in oMin:
                            oMin=str(float(oMin.replace(sKey,EMPTY))*CONVERT[sKey])
                        if sKey in oMax:
                            oMax=str(float(oMax.replace(sKey,EMPTY))*CONVERT[sKey])
                    fMin=float(oMin)
                    fMax=float(oMax)
                    iMean=int((fMax+fMin)/2)
                    dResult[sTarget]=iMean
                    # print(sTarget,iMean)
    return dResult

def WriteOutput(dDict,sPath):
    FILE=open(sPath+"/"+OUTPUT,'w')
    for sKey in sorted(dDict):
        FILE.write("{}{}{}{}".format(sKey,TABULATION,dDict[sKey],ENDLINE))
    FILE.close()

########################################################################
#MAIN
if __name__ == "__main__":
    print("Parsing Viral Zone web pages...")
    dFamily2Length=ParseViralZoneTable()
    print("Write output...")
    WriteOutput(dFamily2Length,sFolder)

########################################################################    
iTime2=time.time()
iDeltaTime=iTime2-iTime1
print("Script done: "+str(iDeltaTime))
