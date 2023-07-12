#! /usr/bin/env python
#>import errorcheck as rcheck
import sys
import def_functions as rfunc

def CheckCompatibilityOfArguements(FalseArguement,TrueArguement):
    if TrueArguement :
        FalseArguement= False
    return FalseArguement

def PrintMolcasLogErrors(filein,lineend):
    linestop = "######.                                           " \
          "                                                  " \
          "                                                  " \
          "                                                  "
    startstrings =["_ERROR_", "not found", "Error", "input error", "Could not find",
                   "not defined", "error", "errorcode",]
    with open(filein, 'r') as fin:
        copy = False
        flag = True
        for line in rfunc.NonBlankLSines(fin):
            matchstring = any(match in line for match in startstrings)
            if matchstring and flag:
                print("!!!!!!!!!!!!!! MOLCAS ENCOUNTERED AN ERROR !!!!!!!!!!!!!!")
                flag=False
                copy = True
            elif linestop in line:
                copy = False
            elif lineend in line:
                copy = False
            elif copy:
                print(line)

    sys.exit()
