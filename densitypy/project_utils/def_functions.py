#! /usr/bin/env python3.6
# Made by Ruben
# >import def_functions as rfunc
import decimal
import json
import sys
from os import path
from subprocess import Popen, PIPE, CalledProcessError

from .logger import logger

def Execute(command):
    # >Executes to command line
    with Popen(command, stdout=PIPE, bufsize=1, universal_newlines=True, shell=True) as p:
        for line in p.stdout:
            logger.info(line, end='')  # process line here
    if p.returncode != 0:
        raise CalledProcessError(p.returncode, p.args)


def execute_command(command: str) -> None:
    """
    Executes a command line instruction and streams the output to a logger.

    :param command: The command to be executed.
    :raises CalledProcessError: If the subprocess exits with a non-zero status.

    Usage:
    >>> execute_command("ls -l")
    """
    logger.info(f"Executing command: {command}")
    # Execute command in a new subprocess
    with Popen(command, stdout=PIPE, bufsize=1, universal_newlines=True, shell=True) as p:
        # Read and log each line of output as it becomes available
        for line in p.stdout:
            logger.info(line.strip())

    # If the subprocess exited with an error, raise an exception
    if p.returncode != 0:
        error_msg = f"Command '{command}' returned non-zero exit status {p.returncode}"
        logger.error(error_msg)
        raise CalledProcessError(p.returncode, p.args, output=p.stdout, stderr=p.stderr)


def ExecuteNoWrite(command):
    # >Executes to command line but does not print
    p = Popen(command, stdout=PIPE, shell=True)
    p_status = p.wait()
    if p_status > 0:
        logger.info("Errors found:: ", p_status)
        sys.exit()


def ExecutePymolcasWithErrorPrint(command, nameofproject):
    # >Executes to command line but does not print
    p = Popen(command, stdout=PIPE, shell=True)
    p_status = p.wait()
    if p_status > 0:
        PrintMolcasLogErrors(nameofproject + ".log", "Timing")


def Execute2(command, goodcall, badcall):
    # >Executes and allows variable prints
    p = Popen(command, stdout=PIPE, shell=True)
    p_status = p.wait()
    if p_status > 0:
        logger.info("Errors found:: ", p_status)
        logger.info(str(badcall))
        sys.exit()
    else:
        logger.info(str(goodcall))


def MakeDirectory(outputdir):
    if path.exists(outputdir):
        execute_command("rm -r " + outputdir)
        execute_command("mkdir " + outputdir)
    else:
        execute_command("mkdir " + outputdir)


def MakeDirectoryNoDelete(outputdir):
    if path.exists(outputdir):
        pass
    else:
        execute_command("mkdir " + outputdir)


def GetValueOfAsString(filename, string, delimiter):
    with open(filename, 'r') as fin:
        for line in fin:
            if string in line:
                option_value = (line.partition(delimiter)[2]).strip()
                option_value = option_value.replace("/", "")
        return option_value


def NonBlankLSines(f):
    for l in f:
        line = l.rstrip()
        if line:
            yield line


def Float_Range(start, stop, step):
    # >Allows non-integer steps
    while start < stop:
        yield format(start, '.1f')  # float(start)
        start += decimal.Decimal(step)


def SaveFile(Filepath, OutPutDir):
    try:
        execute_command("cp " + Filepath + " " + OutPutDir)
    except Exception as e:
        logger.info("Error copying file: " + str(e))


def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def Find(filename, *args):
    directories = [*args]
    foundfile = False
    for searchdirectory in directories:
        if path.exists(searchdirectory + "/" + filename):
            if searchdirectory == ".":
                logger.info(f"Found {filename} inside the current directory")
            else:
                logger.info(f"Found {filename} inside {searchdirectory} directory")
            foundfile = True
            return searchdirectory
    # if not exited by now it means that the file was not found in any of the given directories thus rise error
    if foundfile != True:
        logger.info(f'{filename} not found inside {directories} directories')
        logger.info("exiting...")
        sys.exit()


def GetDipoleValuesAsArray(filename, string, delimiter):
    with open(filename, 'r') as fin:
        value = []
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
    logger.info(filename, extension)
    counter = 1

    while path.exists(filepath):
        filepath = f'{filename}_{counter}{extension}'
        logger.info(counter)
        counter += 1

    return filepath


def remove_json_comments(json_string):
    lines = json_string.split("\n")
    cleaned_lines = [line.split("//")[0] for line in lines]
    cleaned_json = "\n".join(cleaned_lines)
    return cleaned_json


def load_json_file(filename):
    assert filename.endswith(".json") or filename.endswith(".JSON"), "File must be a JSON file"
    with open(filename, 'r') as f:
        cleaned_json = remove_json_comments(f.read())
        return json.loads(cleaned_json)


def CheckCompatibilityOfArguements(FalseArguement, TrueArguement):
    if TrueArguement:
        FalseArguement = False
    return FalseArguement


def PrintMolcasLogErrors(filein, lineend):
    linestop = "######.                                           " \
               "                                                  " \
               "                                                  " \
               "                                                  "
    startstrings = ["_ERROR_", "not found", "Error", "input error", "Could not find",
                    "not defined", "error", "errorcode", ]
    with open(filein, 'r') as fin:
        copy = False
        flag = True
        for line in NonBlankLSines(fin):
            matchstring = any(match in line for match in startstrings)
            if matchstring and flag:
                print("!!!!!!!!!!!!!! MOLCAS ENCOUNTERED AN ERROR !!!!!!!!!!!!!!")
                flag = False
                copy = True
            elif linestop in line:
                copy = False
            elif lineend in line:
                copy = False
            elif copy:
                print(line)

    sys.exit()
