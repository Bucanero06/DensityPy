#! /usr/bin/env python3.6
#Made by Ruben
#>import def_functions as rfunc

# import subprocess #commented to test!!!!!!!!!!!!!!

from subprocess import Popen, PIPE, CalledProcessError
import decimal
from os import path, system
import sys
import pandas as pd
from errorcheck import PrintMolcasLogErrors


def Execute(command):
    #>Executes to command line
    with Popen(command, stdout=PIPE, bufsize=1, universal_newlines=True, shell=True) as p:
        for line in p.stdout:
            print(line, end='')  # process line here
    if p.returncode != 0:
        raise CalledProcessError(p.returncode, p.args)

	
def ExecuteNoWrite(command):
    #>Executes to command line but does not print
    p = Popen(command, stdout=PIPE, shell=True)
    p_status = p.wait()
    if p_status > 0:
        print("Errors found:: ", p_status)
        sys.exit()

def ExecutePymolcasWithErrorPrint(command, nameofproject):
    #>Executes to command line but does not print
    p = Popen(command, stdout=PIPE, shell=True)
    p_status = p.wait()
    if p_status > 0:
        PrintMolcasLogErrors(nameofproject + ".log", "Timing")
        
def Execute2(command,goodcall,badcall):
    # >Executes and allows variable prints
    p = Popen(command, stdout=PIPE, shell=True)
    p_status = p.wait()
    if p_status > 0:
        print("Errors found:: ", p_status)
        print(str(badcall))
        sys.exit()
    else:
        print(str(goodcall))

def MakeDirectory(outputdir):
    if path.exists(outputdir):
        Execute("rm -r " + outputdir)
        Execute("mkdir " + outputdir)
    else:
        Execute("mkdir " + outputdir)
def MakeDirectoryNoDelete(outputdir):
    if path.exists(outputdir):
        pass
    else:
        Execute("mkdir " + outputdir)

def GetValueOfAsString(filename,string,delimiter):
    with open(filename, 'r') as fin:
        for line in fin:
            if string in line:
                option_value = (line.partition(delimiter)[2]).strip()
                option_value=option_value.replace("/","")
        return option_value

def NonBlankLSines(f):
    for l in f:
        line = l.rstrip()
        if line:
            yield line

def Float_Range(start, stop, step):
    #>Allows non-integer steps
    while start < stop:
        yield  format(start, '.1f') #float(start)
        start += decimal.Decimal(step)

def SaveFile(FileName, OutPutDir):
    Execute("cp " + FileName + " " + OutPutDir )

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def Find(filename, *args):
    directories=[*args]
    foundfile = False
    for searchdirectory in directories:
        if path.exists(searchdirectory + "/" + filename):
            if searchdirectory == ".":
                print("Found " + str(filename) + " inside the current directory")
            else:
                print("Found " + str(filename) +  " inside " + str(searchdirectory) + " directory")
            foundfile = True
            return searchdirectory
            exit()
    # if not exited by now it means that the file was not found in any of the given directories thus rise error
    if foundfile != True:
        print(str(filename) + " not found inside " + str(directories) + "\n exiting...")
        sys.exit()

def GetDipoleValuesAsArray(filename,string,delimiter):
    with open(filename, 'r') as fin:
        value=[]
        for line in fin:
            if string in line:
                option_value = (line.partition(delimiter)[2]).strip()
                value.append(option_value)
        return value

def File_Lenth(filename):
    with open(filename) as fin:
        for i, l in enumerate(fin):
            pass
    return i + 1

def uniquify(filepath):
    filename, extension = path.splitext(filepath)
    print(filename , extension)
    counter = 1

    while path.exists(filepath):
        filepath = f'{filename}_{counter}{extension}'
        print(counter)
        counter += 1

    return filepath