#! /usr/bin/env python3.10
import decimal
import json
import os
import re
import shutil
import sys
from contextlib import contextmanager
from os import path
from subprocess import Popen, PIPE, CalledProcessError

import h5py as h5py
import numpy as np

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
                    _logger.info(f' ---> {line}')
            for line in p.stderr:
                line = line.strip()
                stderr_output_list.append(line)
                if write_stream:
                    _logger.error(f' ---> {line}')
        return_code = p.returncode

        stdout_output_list, stderr_output_list = json.dumps(stdout_output_list), json.dumps(stderr_output_list)
        if return_code != 0:
            raise CalledProcessError(return_code, p.args, output=stdout_output_list,
                                     stderr=stderr_output_list)

        return return_code, json.loads(stdout_output_list), json.loads(stderr_output_list)

    except Exception as e:
        raise e


@contextmanager
def change_directory(new_directory):
    """Context manager to change the working directory temporarily."""
    current_directory = os.getcwd()
    os.chdir(new_directory)
    try:
        yield
    finally:
        os.chdir(current_directory)


def extract_datasets_from_h5_to_csv(h5_filepath, dataset_mapping):
    """
    Extract specific datasets from an HDF5 file and save them to CSV files.

    Parameters:
    - h5_filepath: Path to the input HDF5 file.
    - dataset_mapping: Dictionary where keys are dataset names in the HDF5 file
                       and values are the output filenames for CSV files.

    Note:
        - it will flatten the data in the dataset to a 2D array if the dataset is a 3D array
            - 3D arrays of shape (x, y, z) → 2D arrays of shape (x*y, z)
            - 4D arrays of shape (w, x, y, z) → 2D arrays of shape (wxy, z)
            - ... and so on.
    """
    with h5py.File(h5_filepath, 'r') as h5_file:
        for dataset_name, output_file in dataset_mapping.items():
            try:
                if dataset_name in h5_file:
                    logger.debug(f"Extracting dataset {dataset_name} from {h5_filepath} to {output_file}")
                    data = flatten_ND_to_2D(h5_file[dataset_name][:]) # Flatten the data to 2D
                    # Save the data as a CSV file
                    np.savetxt(output_file, data, delimiter=',')
                else:
                    logger.warning(f"Dataset {dataset_name} not found in {h5_filepath}")
            except Exception as e:
                logger.error(f"Error extracting dataset {dataset_name}: {e}")


def flatten_ND_to_2D(array):
    """Flatten a multi-dimensional array keeping the last dimension and format its values."""
    if array.ndim > 1:
        flat = array.reshape(-1, array.shape[-1])  # Flatten all dimensions except the last one
    else:
        flat = array  # If already 1D, no reshaping needed
    return flat


def remove_spaces_and_to_lowercase(key: str) -> str:
    """
    Normalize a string key by converting all characters to lowercase
    and removing all spaces and underscores.

    :param key: The string key to be normalized.
    :type key: str
    :return: The normalized key.
    :rtype: str

    Usage::

        normalized_key = normalize_key('Hello_World')
        print(normalized_key)  # Output: 'helloworld'
    """
    translation_table = str.maketrans('', '', ' _')
    return key.lower().translate(translation_table)


def recursively_normalize_dict_keys(d: dict) -> dict:
    """
    Normalize all keys in a dictionary. If a value is a dictionary,
    normalize its keys recursively.

    :param d: The dictionary whose keys are to be normalized.
    :type d: dict
    :return: A new dictionary with normalized keys.
    :rtype: dict

    Usage::

        data = {
            'Hello World': 1,
            'Good_Day': {
                'Inner_Key': 2
            }
        }
        normalized_data = normalize_dict(data)
        print(normalized_data)  # Output: {'helloworld': 1, 'goodday': {'innerkey': 2}}
    """
    return {
        remove_spaces_and_to_lowercase(key): recursively_normalize_dict_keys(value) if isinstance(value,
                                                                                                  dict) else value
        for key, value in d.items()
    }


def validate_and_correct_dictionary(dict_to_verify: dict, default_config: dict, log=None) -> tuple:
    """
    Verifies that the given JSON config contains all the required keys as per the default config.
    If not, adds the missing keys with the default values.
    Returns a log of warnings and errors in case of missing or unknown keys.

    :param dict_to_verify: The JSON configuration to verify.
    :type dict_to_verify: dict
    :param default_config: The default configuration used for verification.
    :type default_config: dict
    :param log: A list of log messages. Default is an empty list.
    :type log: list, optional
    :return: A tuple of the verified JSON configuration and the log of warnings/errors.
    :rtype: tuple

    Usage::

        json_config = {
            "key1": "value1",
            "key2": {
                "subkey1": "subvalue1"
            }
        }

        default_config = {
            "key1": "value1",
            "key2": {
                "subkey1": "subvalue1",
                "subkey2": "subvalue2"
            },
            "key3": "value3"
        }

        verified_config, log = verify_configuration_keys(json_config, default_config)
        # verified_config will contain all keys from default_config, and log will contain warning/error messages
    """
    if log is None:
        log = []
    if not dict_to_verify:
        dict_to_verify = {}

    for key, value in default_config.items():
        if key not in dict_to_verify:
            dict_to_verify[key] = value
            log.append(f'Warning: Missing key "{key}" has been added with default value "{value}"')
        elif isinstance(value, dict):
            if isinstance(dict_to_verify[key], dict):
                dict_to_verify[key], log = validate_and_correct_dictionary(dict_to_verify[key], value, log)
            else:
                log.append(f'Error: Key "{key}" should contain a dictionary but found "{type(dict_to_verify[key])}"')
        else:
            if not isinstance(dict_to_verify[key], type(value)):
                log.append(
                    f'Error: Key "{key}" should be of type "{type(value)}" but found "{type(dict_to_verify[key])}"')

    for key in dict_to_verify.keys():
        if key not in default_config:
            log.append(f'Error: Unknown key "{key}" found in config')

    return dict_to_verify, log


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
        logger.info(f"Copying file: {input} to {output}")
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


def delete_files_or_directories(*paths):
    """
    Deletes the files or directories specified by the given paths.

    :param paths: The paths of files or directories to be deleted.
    :type paths: str
    :raises FileNotFoundError: If the given path does not exist.

    Usage::

        delete_files_or_directories('/path/to/file', '/path/to/directory')  # Deletes specified file and directory
    """
    for paths in paths:
        if os.path.isfile(paths):
            os.remove(paths)  # Delete the file
        elif os.path.isdir(paths):
            shutil.rmtree(paths)  # Delete the directory and its contents
        else:
            raise FileNotFoundError(f"Path '{paths}' does not exist.")


def make_directory(output_dir: str, delete_if_exists: bool = False):
    """
    Creates a directory at the specified path. If the directory already exists and 'delete_if_exists' is True,
    it deletes the existing directory before creating a new one.

    :param output_dir: The path where the directory is to be created.
    :type output_dir: str
    :param delete_if_exists: If True, deletes the directory at 'output_dir' if it exists. Default is False.
    :type delete_if_exists: bool

    Usage::

        make_directory('/path/to/directory', delete_if_exists=True)  # Creates directory, deleting existing one if necessary
    """
    if delete_if_exists and os.path.exists(output_dir):  delete_files_or_directories(output_dir)
    os.makedirs(output_dir, exist_ok=True)  # Create the directory


class change_directory_manager:
    def __init__(self, new_dir):
        self.new_dir = new_dir
        self.original_dir = os.getcwd()

    def __enter__(self):
        os.chdir(self.new_dir)

    def __exit__(self, type, value, traceback):
        os.chdir(self.original_dir)
