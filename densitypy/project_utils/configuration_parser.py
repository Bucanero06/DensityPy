#! /usr/bin/env python3.10
import json
from os import path

from .def_functions import load_json_file
from ..Default_Settings.default_config import DEFAULT_CONFIG_FILE_PATH, DEFAULT_CONFIG_CONTENT

from densitypy.project_utils.logger import setup_logger

logger = setup_logger(__name__.split('.')[-1])


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
    translation_table = str.maketrans('', '', ' _')
    return key.lower().translate(translation_table)


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
    return {
        normalize_key(key): normalize_dict(value) if isinstance(value, dict) else value
        for key, value in d.items()
    }

import json


def verify_configuration_keys(json_config: dict, default_config: dict, log: list = []) -> tuple:
    """
    Verifies that the given JSON config contains all the required keys as per the default config.
    If not, adds the missing keys with the default values.
    Returns a log of warnings and errors in case of missing or unknown keys.

    :param json_config: The JSON configuration to verify.
    :type json_config: dict
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
    if not json_config:
        json_config = {}

    for key, value in default_config.items():
        if key not in json_config:
            json_config[key] = value
            log.append(f'Warning: Missing key "{key}" has been added with default value "{value}"')
        elif isinstance(value, dict):
            if isinstance(json_config[key], dict):
                json_config[key], log = verify_configuration_keys(json_config[key], value, log)
            else:
                log.append(f'Error: Key "{key}" should contain a dictionary but found "{type(json_config[key])}"')
        else:
            if not isinstance(json_config[key], type(value)):
                log.append(f'Error: Key "{key}" should be of type "{type(value)}" but found "{type(json_config[key])}"')

    for key in json_config.keys():
        if key not in default_config:
            log.append(f'Error: Unknown key "{key}" found in config')

    return json_config, log


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

    default_config = normalize_dict(DEFAULT_CONFIG_CONTENT)
    config, log = verify_configuration_keys(config, default_config)

    for message in log:
        logger.warning(message)


    assert config, f'Could not load configuration file {config_file_path}'
    from pprint import pprint
    pprint(config)
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
