#! /usr/bin/env python3.10

from densitypy.project_utils.logger import setup_logger

logger = setup_logger(__name__.split('.')[-1])


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
