#! /usr/bin/env python3.10
import decimal
import json
import re
import sys
from os import path
from subprocess import Popen, PIPE, CalledProcessError

from .logger import logger


def execute_command(command: str, write_stream=True) -> None:
    """
    Executes a command line instruction and streams the output to a logger.

    :param command: The command to be executed.
    :param write_stream: A boolean flag indicating if the output should be streamed.
    :type command: str
    :type write_stream: bool, optional
    :raises CalledProcessError: If the subprocess exits with a non-zero status.

    Usage::

    >>> execute_command("ls -l")
    """
    logger.info(f"Executing command: {command}")
    # Execute command in a new subprocess
    with Popen(command, stdout=PIPE, bufsize=1, universal_newlines=True, shell=True) as p:
        # Read and log each line of output as it becomes available
        # todo accept errors and warning
        if write_stream:
            for line in p.stdout:
                logger.info(line.strip())

    # If the subprocess exited with an error, raise an exception
    if p.returncode != 0:
        error_msg = f"Command '{command}' returned non-zero exit status {p.returncode}"
        logger.error(error_msg)
        raise CalledProcessError(p.returncode, p.args, output=p.stdout, stderr=p.stderr)


def execute_pymolcas_with_error_print(command, nameofproject):
    """
    Executes a command line instruction but does not print the output.
    If the command execution fails, it prints the errors in the log file.

    :param command: The command to be executed.
    :param nameofproject: The name of the project.
    :type command: str
    :type nameofproject: str
    """
    p = Popen(command, stdout=PIPE, shell=True)
    p_status = p.wait()
    if p_status > 0:
        print_molcas_log_errors(nameofproject + ".log", "Timing")


def make_directory(outputdir):
    """
    Creates a new directory. If the directory already exists, it removes it first.

    :param outputdir: The path to the directory to be created.
    :type outputdir: str
    """
    if path.exists(outputdir):
        execute_command(f'rm -r {outputdir}')
        execute_command(f'mkdir {outputdir}')
    else:
        execute_command(f'mkdir {outputdir}')


def make_directory_no_delete(outputdir):
    """
    Creates a new directory. If the directory already exists, it does nothing.

    :param outputdir: The path to the directory to be created.
    :type outputdir: str
    """
    if path.exists(outputdir):
        pass
    else:
        execute_command("mkdir " + outputdir)


def get_value_of_as_string(filename, string, delimiter):
    """
    Reads a file and returns the value of a string in the file.

    :param filename: The path to the file.
    :param string: The string to search for.
    :param delimiter: The delimiter that separates the string and its value.
    :type filename: str
    :type string: str
    :type delimiter: str
    :return: The value of the string in the file.
    :rtype: str
    """
    with open(filename, 'r') as fin:
        for line in fin:
            if string in line:
                option_value = (line.partition(delimiter)[2]).strip()
                option_value = option_value.replace("/", "")
        return option_value


def non_blank_l_sines(f):
    """
    A generator function that yields non-blank lines from a file.

    :param f: The file to read.
    :type f: file
    :yields: Each non-blank line in the file.
    :rtype: str
    """
    for l in f:
        line = l.rstrip()
        if line:
            yield line


def float_range(start, stop, step):
    """
    A generator function that yields a range of floating point numbers.

    :param start: The start of the range.
    :param stop: The end of the range.
    :param step: The step size.
    :type start: float
    :type stop: float
    :type step: float
    :yields: Each value in the range.
    :rtype: float
    """
    while start < stop:
        yield format(start, '.1f')  # float(start)
        start += decimal.Decimal(step)


def copy_file_to(input, output):
    """
    Copies a file to a new location.

    :param input: The path of the file to be copied.
    :param output: The path to copy the file to.
    :type input: str
    :type output: str
    """
    try:
        execute_command(f'cp {input} {output}')
    except Exception as e:
        logger.info("Error copying file: " + str(e))


def file_len(file_path):
    """
    Returns the number of lines in a file.

    :param file_path: The path to the file.
    :type file_path: str
    :return: The number of lines in the file.
    :rtype: int
    """
    with open(file_path) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def find(filename, *args):
    """
    Searches for a file in multiple directories.

    :param filename: The name of the file to search for.
    :param args: The directories to search in.
    :type filename: str
    :type args: str
    :return: The first directory where the file was found.
    :rtype: str
    """
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


def get_dipole_values_as_array(filename, string, delimiter):
    """
    Reads a file and returns the values of a string in the file as an array.

    :param filename: The path to the file.
    :param string: The string to search for.
    :param delimiter: The delimiter that separates the string and its values.
    :type filename: str
    :type string: str
    :type delimiter: str
    :return: The values of the string in the file.
    :rtype: list
    """
    with open(filename, 'r') as fin:
        value = []
        for line in fin:
            if string in line:
                option_value = (line.partition(delimiter)[2]).strip()
                value.append(option_value)
        return value


def file_lenth(filename):
    """
    Returns the number of lines in a file.

    :param filename: The path to the file.
    :type filename: str
    :return: The number of lines in the file.
    :rtype: int
    """
    with open(filename) as fin:
        for i, l in enumerate(fin):
            pass
    return i + 1


def uniquify(filepath):
    """
    Makes a file path unique by appending a counter to the file name.

    :param filepath: The initial file path.
    :type filepath: str
    :return: A unique file path.
    :rtype: str
    """
    filename, extension = path.splitext(filepath)
    logger.info(filename, extension)
    counter = 1

    while path.exists(filepath):
        filepath = f'{filename}_{counter}{extension}'
        logger.info(counter)
        counter += 1

    return filepath


def remove_json_comments(json_string):
    """
    Removes the comments from a JSON string.

    :param json_string: The JSON string.
    :type json_string: str
    :return: The cleaned JSON string.
    :rtype: str
    """
    lines = json_string.split("\n")
    cleaned_lines = [line.split("//")[0] for line in lines]
    cleaned_json = "\n".join(cleaned_lines)
    return cleaned_json


def load_json_file(filename):
    """
    Loads a JSON file.

    :param filename: The path to the JSON file.
    :type filename: str
    :return: The content of the JSON file.
    :rtype: dict
    """
    assert filename.endswith(".json") or filename.endswith(".JSON"), "File must be a JSON file"
    with open(filename, 'r') as f:
        cleaned_json = remove_json_comments(f.read())
        return json.loads(cleaned_json)


def check_compatibility_of_arguements(FalseArguement, TrueArguement):
    """
    Checks the compatibility of two arguments.

    :param FalseArguement: The false argument.
    :param TrueArguement: The true argument.
    :type FalseArguement: bool
    :type TrueArguement: bool
    :return: The compatibility status.
    :rtype: bool
    """
    if TrueArguement:
        FalseArguement = False
    return FalseArguement


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


def natural_sort(iterable, key=None, reverse=False):
    """
    Sorts the given iterable in a natural order.

    This function is a key-function to the built-in `sorted` function and can be
    used as a drop-in replacement for it.

    A natural sort, also known as an alphanumeric sort, is a sorting method that orders strings containing numbers in
    a way that considers the numerical value of the digits rather than treating the entire string as a sequence of
    characters. In other words, it sorts strings with numbers in a way that reflects their numerical order.

    :param iterable: The iterable to be sorted.
    :param key: A callable used to extract a comparison key from each element in the iterable.
    :param reverse: If set to True, the iterable will be sorted in descending order.
    :type iterable: iterable
    :type key: callable, optional
    :type reverse: bool, optional
    :return: A new list containing the sorted elements from the iterable.
    :rtype: list

    Usage::
        >>> natural_sort(['2 ft', '10 ft', '1 ft'])
        ['1 ft', '2 ft', '10 ft']
    """

    def __float_convert(match):
        try:
            return float(match.group())
        except ValueError:
            return match.group()

    if key is None:
        key = lambda x: x
    else:
        key = lambda x: (__float_convert(match) for match in re.finditer(r'\d+|\D+', key(x)))

    return sorted(iterable, key=key, reverse=reverse)
