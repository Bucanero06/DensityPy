#! /usr/bin/env python3.10
import json
from os import path

from .def_functions import load_json_file
from .logger import logger
from ..Default_Settings.default_config import DEFAULT_CONFIG_FILE_PATH, DEFAULT_CONFIG_CONTENT


def normalize_key(key: str) -> str:
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
    return key.lower().replace(' ', '').replace('_', '')


def normalize_dict(d: dict) -> dict:
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
    new_dict = {}
    for key, value in d.items():
        if isinstance(value, dict):
            # If the value is a dictionary, normalize its keys recursively
            new_dict[normalize_key(key)] = normalize_dict(value)
        else:
            # If the value is not a dictionary, just copy the value
            new_dict[normalize_key(key)] = value
    return new_dict


def parse_configuration_file(config_file_path=None):
    """
    Parses a configuration file and returns a dictionary with the configuration.

    If a configuration file path is provided and the file exists, it will be used to load the configuration.
    If the provided file path does not exist, or if no path is provided, a default configuration file will be created
    and a warning message will be logged.

    :param config_file_path: The path to the configuration file. Default is None.
    :type config_file_path: str, optional
    :return: A dictionary representing the configuration.
    :rtype: dict
    :raises AssertionError: If the configuration file cannot be loaded.

    Usage::

        config = parse_configuration_file('config.json')  # Load configuration from 'config.json'
    """
    load_user_configuration = config_file_path and path.isfile(config_file_path)

    if load_user_configuration:
        logger.info(f'Found {config_file_path}')
    else:
        logger.warning(
            f'{config_file_path if config_file_path else "No configuration file"} found. '
            f'Making default configuration file with name {DEFAULT_CONFIG_FILE_PATH}. '
            'Check/Edit configuration file before running.'
        )
        write_example_configuration_file(DEFAULT_CONFIG_FILE_PATH)
    #
    config = load_json_file(config_file_path if load_user_configuration else DEFAULT_CONFIG_FILE_PATH)
    config = normalize_dict(config)
    assert config, f'Could not load configuration file {config_file_path}'
    return config


def write_example_configuration_file(file_name: str = DEFAULT_CONFIG_FILE_PATH):
    """
    Writes an example configuration file.

    This function creates a configuration file with default settings in the current working directory.
    If a file name is not provided, it uses a default name.

    :param file_name: The name of the file to be created. Default is DEFAULT_CONFIG_FILE_PATH.
    :type file_name: str, optional
    :raises Exception: Logs an error if the example configuration file cannot be written.

    Usage::

        write_example_configuration_file('example_config.json') # Creates 'example_config.json' with default settings
    """

    try:
        # Write Configuration FIle with Defaults
        with open(file_name, 'w') as configfile:
            json.dump(DEFAULT_CONFIG_CONTENT, configfile, indent=4)
        logger.info(f'Wrote example configuration file {file_name}')
    except Exception as e:
        logger.error(f'Could not write example configuration file {file_name}')
        logger.error(e)
