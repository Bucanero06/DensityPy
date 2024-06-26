#! /usr/bin/env python3.10
import json
import sys
from subprocess import Popen, PIPE, CalledProcessError

from densitypy.project_utils.file_parsing import non_blank_l_sines
from densitypy.project_utils.logger import setup_logger

logger = setup_logger(__name__.split('.')[-1])


# todo i want to be able to pass a re output handler for user defined parsing
def execute_command(command: str, write_stream=True, _logger=None) -> (int, list, list):
    """
    Execute a shell command and capture the standard output and standard error streams.

    :param command: The shell command to execute.
    :type command: str
    :param write_stream: Flag indicating whether to write the output streams to the log, defaults to True.
    :type write_stream: bool, optional
    :param _logger: A logger instance to use for logging, defaults to None, in which case the global `logger` is used.
    :type _logger: logging.Logger, optional
    :return: The return code of the command, and the standard output and standard error captured as lists.
    :rtype: tuple(int, list, list)
    :raises CalledProcessError: If the command returns a non-zero exit status.
    """
    if _logger is None:
        _logger = logger

    _logger.info(f"Executing command: {command}")
    stdout_output_list = []
    stderr_output_list = []
    try:
        with Popen(command, stdout=PIPE, stderr=PIPE, bufsize=1, universal_newlines=True, shell=True) as p:
            for line in p.stdout:
                line = line.strip()
                stdout_output_list.append(line)
                if write_stream:
                    _logger.info(line)
            for line in p.stderr:
                line = line.strip()
                stderr_output_list.append(line)
                if write_stream:
                    _logger.error(line)
        return_code = p.returncode

        stdout_output_list, stderr_output_list = json.dumps(stdout_output_list), json.dumps(stderr_output_list)
        if return_code != 0:
            raise CalledProcessError(return_code, p.args, output=stdout_output_list,
                                     stderr=stderr_output_list)

        return return_code, json.loads(stdout_output_list), json.loads(stderr_output_list)

    except Exception as e:
        raise e




def print_molcas_log_errors(filein, lineend):
    """
    Prints the errors in a Molcas log file.

    :param filein: The input log file.
    :param lineend: The line ending marker.
    :type filein: str
    :type lineend: str
    """
    linestop = "######.                                           " \
               "                                                  " \
               "                                                  " \
               "                                                  "
    startstrings = ["_ERROR_", "not found", "Error", "input error", "Could not find",
                    "not defined", "error", "errorcode", ]
    with open(filein, 'r') as fin:
        copy = False
        flag = True
        for line in non_blank_l_sines(fin):
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
