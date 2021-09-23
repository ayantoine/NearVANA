# coding: utf-8
"""Python3.6"""

import time
from optparse import OptionParser
import pandas
import requests
import os

sCurrentVersionScript="v1"
iTime1=time.time()
########################################################################
'''
V1-2021/09/23
Download and convert the VMR files directly from the ICTVwebsite

python GetICTV-VMRfile.py -f FOLDER
FOLDER: output folder where to store the table
'''
########################################################################
#CONSTANT
OUTPUT_SUFFIX=".tsv"

VMR_EXPLICIT="Virus Metadata Repository"
DOWNLOAD="download"
HREF="href"

ICTV_VMR_MAINPAGE="https://talk.ictvonline.org/taxonomy/vmr/"
ICTV_BASENAME="https://talk.ictvonline.org/"

ICTV_VMR_XLSX="ITCV-VMR.xlsx"

ENDLINE="\n"
TABULATION="\t"
DOUBLEQUOTE="\""


########################################################################
#Options
parser = OptionParser()
parser.add_option("-f","--folder", dest="folder")

(options, args) = parser.parse_args()

sFolder=options.folder
if not sFolder:
    exit("Error : no folder -f defined, process broken")

########################################################################
#Function
def DownloadFile(sFolder):
    oMainPage=requests.get(ICTV_VMR_MAINPAGE)
    for sNewLine in oMainPage.text.split(ENDLINE):
        if VMR_EXPLICIT in sNewLine:
            tLine=sNewLine.split(DOUBLEQUOTE)
            sLink=tLine[1]
            break
    print("Latest VMR available at : "+sLink)
    oSubPage=requests.get(sLink)
    for sNewLine in oSubPage.text.split(ENDLINE):
        if DOWNLOAD in sNewLine and HREF in sNewLine:
            tLine=sNewLine.split(DOUBLEQUOTE)
            sLink=ICTV_BASENAME+tLine[1]
            break
    print("Download file : "+sLink)
    oRequest=requests.get(sLink,allow_redirects=True)
    open(sFolder+"/"+ICTV_VMR_XLSX,"wb").write(oRequest.content)

def ConvertFile(sFolder):
    oXlsx=pandas.ExcelFile(ICTV_VMR_XLSX)
    sFirstSheet=oXlsx.sheet_names[0]
    oXlsx=pandas.read_excel(ICTV_VMR_XLSX, sFirstSheet, index_col=None)
    oXlsx.to_csv(sFirstSheet+OUTPUT_SUFFIX,sep="\t",encoding="utf-8",index=False,line_terminator=ENDLINE)

########################################################################
#MAIN
if __name__ == "__main__":
    print("Search and Download VMC on ICTV...")
    DownloadFile(sFolder)
    print("Convert file into tsv format...")
    ConvertFile(sFolder)
    print("Remove downloaded file...")
    os.remove(sFolder+"/"+ICTV_VMR_XLSX)

########################################################################    
iTime2=time.time()
iDeltaTime=iTime2-iTime1
print("Script done: "+str(iDeltaTime))
