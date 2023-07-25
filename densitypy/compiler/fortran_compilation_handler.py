"""This code is defined in Python, specifically intended to facilitate the process of compiling Fortran code. It
provides logging facilities to capture any errors or warnings that may occur during the compilation process. It
contains two key methods:

1. `process_compiler_output(output: str)`: This function is used to parse the output of the Fortran compiler to
identify errors and warnings, and log them in a more readable format. It does this by looping over each line of the
output, using a regular expression to extract relevant information, and then storing the errors and warnings in
defaultdict objects.

2. `compile_fortran_code(directory: str, make_flags: str = '')`: This function is used to compile the Fortran code in
the given directory and handle any errors or warnings. It changes to the given directory, compiles the code,
and then changes back to the original directory. If any errors are encountered during the process, they are logged
and further processed by `process_compiler_output(output: str)` function.

This code also imports various modules for usage:

- `os`: To interact with the OS, specifically for changing directories and getting current working directory. - `re`:
For regular expressions. - `defaultdict`: A type of Python's built-in dictionary that allows you to specify a default
value for keys that have not yet been assigned a value. - `CalledProcessError`: A specific type of Python error,
raised when a process run by the subprocess module (or by the command line) returns a non-zero exit status. -
`execute_command`: A custom function from a local module to execute a shell command. - `setup_logger`: A custom
function from a local module to set up a logger for capturing and handling log messages.

The imported `setup_logger` function is used to create a logger, which is then used within these functions to log any
relevant information, such as successful compilation, or any errors or warnings that may have occurred.

Note: This is a fairly advanced script, and it assumes a certain file structure and that certain utilities (like a
'make' command and a proper setup Makefile) are available in the directory. The script also assumes the presence of
certain locally defined modules (`densitypy.project_utils.def_functions.execute_command` and
`densitypy.project_utils.logger.setup_logger`).

"""

import os
import re
from collections import defaultdict
from subprocess import CalledProcessError

import pandas as pd

from densitypy.project_utils.def_functions import execute_command_with_capture
from densitypy.project_utils.logger import setup_logger

logger = setup_logger("compiler")


def print_in_human_readable_format(remarks, warnings, errors, others):
    if remarks:
        for file_path, remark_list in remarks.items():
            logger.remarks(f"File: {file_path}")
            for line_number, message_type, code, message, variable, code_reference in remark_list:
                logger.remarks(f"Line {line_number}{code}:[{variable}] {message} ")

    if warnings:
        for file_path, warning_list in warnings.items():
            logger.warning(f"File: {file_path}")
            for line_number, message_type, code, message, variable, code_reference in warning_list:
                logger.warning(f"Line {line_number}{code}:[{variable}] {message} ")

    if others:
        for file_path, others_list in others.items():
            for line_number, message_type, code, message, variable, code_reference in others_list:
                logger.unknown(f"{message}")

    if errors:
        for file_path, error_list in errors.items():
            logger.error(f"File: {file_path}")
            for line_number, message_type, code, message, variable, code_reference in error_list:
                logger.error(f"Line {line_number}{code}:[{variable}] {message} ")



def process_compiler_output(output: str, return_levels: bool = False):
    # Dictionaries to store errors, warnings, remarks, and others
    errors = defaultdict(list)
    warnings = defaultdict(list)
    remarks = defaultdict(list)
    others = defaultdict(list)

    # Create an empty list to store the dictionaries
    data_list = []

    # regex pattern to match error, warning, and remark messages
    pattern = r"^(.+)\((\d+)\):\s+(.+?)\s+#(\d+):\s+(.+?)\.\s+\[(.+?)\]"
    DUMMY_LINE_COUNTER = None

    output_list = output.split("\COMMA")
    for index, line in enumerate(output_list):
        match = re.match(pattern, line)
        if match:
            DUMMY_LINE_COUNTER = 0

            file_path, line_number, message_type, code, message, variable = match.groups()
            line_number = int(line_number)
            code = int(code)
            code_reference = output_list[index + 1]

            # Save match into appropriate category
            if "warning" in message_type.lower():
                warnings[file_path].append((line_number, message_type, code, message, variable, code_reference))
            elif "error" in message_type.lower():
                errors[file_path].append((line_number, message_type, code, message, variable, code_reference))
            elif "remark" in message_type.lower():
                remarks[file_path].append((line_number, message_type, code, message, variable, code_reference))
        else:
            if DUMMY_LINE_COUNTER == 0:
                DUMMY_LINE_COUNTER = 1
            elif DUMMY_LINE_COUNTER == 1:
                pass
            else:
                DUMMY_LINE_COUNTER = None
                # save lines which don't match the pattern in others
                file_path, line_number, message_type, code, message, variable, code_reference = 'Unknown', pd.NA, 'other', pd.NA, line, 'Unknown', 'Unknown'
                others[file_path].append((line_number, message_type, code, message, variable, code_reference))

        if DUMMY_LINE_COUNTER != 1:
            # Append the extracted data as a dictionary to the list
            data_list.append({
                'File': file_path,
                'Line_Number': line_number,
                'Type': message_type,
                'Code': code,
                'Message': message,
                'Variable': variable,
                'Code_Reference': code_reference
            })
    # Create the DataFrame from the list of dictionaries
    df = pd.DataFrame(data_list)

    # Sort the DataFrame by line number (with -1 values moved to the bottom)
    df = df.sort_values(by=['Line_Number'], ascending=True, na_position='last').astype(
        {'Line_Number': 'Int64', 'Code': 'Int64'})
    df = df.reset_index(drop=True)

    # Print output in human readable format
    print_in_human_readable_format(remarks=remarks, warnings=warnings, errors=errors, others=others)

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # do not restrict width of columns
    pd.set_option('display.max_colwidth', None)
    # Print the DataFrame
    print(df[["Variable", "Type", "Message", "Code_Reference"]])
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    # Return appropriate output based on flag
    if return_levels:
        return errors, warnings, remarks, others, df
    else:
        return df


def compile_fortran_code(directory: str, make_flags: str = ''):
    """
    Compile the Fortran code in the given directory and handle any errors or warnings.

    :param directory: The directory containing the Fortran code to be compiled.
    :param make_flags: Additional flags to be passed to the Makefile.
    :type directory: str.
    :type make_flags: str, optional.
    """
    # Use the given directory and get the exact path needed to not cause any errors
    directory = os.path.abspath(directory)
    # Compile the code and handle any potential errors
    try:
        # Compile the code and capture output and errors
        execute_command_with_capture(f"make -C {directory} {make_flags}", write_stream=False)
        logger.info(f"Successfully compiled or already found compiled code in {directory}")
    except CalledProcessError as e:
        logger.error(f"Failed to compile code in {directory}")
        # Read the content of the error stream and pass it as a string to process_compiler_output()
        process_compiler_output(e.stderr)

    exit()  # fixme Exit the program for now as im developing

# from pydantic import BaseModel, Field, validator
# import os
#
# class FortranCompilerConfig(BaseModel):
#     directory: str = Field(..., description="The directory containing the Fortran code to be compiled.")
#     make_flags: str = Field(default='', description="Additional flags to be passed to the Makefile.")
#
#     @validator('directory')
#     def directory_exists(cls, v):
#         if not os.path.isdir(v):
#             raise ValueError(f"The directory {v} does not exist.")
#         return v
#
#     @validator('make_flags')
#     def validate_make_flags(cls, v):
#         if not isinstance(v, str):
#             raise ValueError("make_flags must be a string.")
#         return v
#
#     class Config:
#         validate_assignment = True
# def compile_fortran_code_with_config(config: FortranCompilerConfig):
#     """
#     Compile the Fortran code using a FortranCompilerConfig object.
#
#     :param config: A FortranCompilerConfig object containing directory and make_flags.
#     """
#     # Use the config object for compilation
#     compile_fortran_code(config.directory, config.make_flags)
#
