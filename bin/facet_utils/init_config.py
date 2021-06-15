#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Alex Stewart"

import json
import os.path
from pathlib import Path


parent_dir = "%s/" % Path(os.path.realpath(__file__)).parent
config_filename = "PROG_OUTDIR_CONFIG"


def validate_config_load(config_dict):
    if try_config_key(config_dict, "db_name") and try_config_key(config_dict, "outdir_name"):
        return config_dict
    else:
        rewrite = input("One of the required config entries doesn't exist!\n"
                        "Would you like to rewrite a default config file? (y/n) ")
        if rewrite == 'y':
            write_default_config()
            return get_config()

        else:
            print("Please fix the config file!")
            print("Exiting....")
            exit()


def get_config():
    try:
        with open(parent_dir+config_filename, 'r') as config:
            return json.load(config)

    except Exception as e:
        print("There was a problem loading the %s file in %s:" % (config_filename, parent_dir))
        print(e, "\n")
        rewrite = input("Would you like to rewrite a default config file? (y/n) ")
        if rewrite == 'y':
            write_default_config()
            return get_config()
        else:
            print("\nPlease check this file to make sure it exists, you can access it, and its contents are okay\n")
            print("Exiting....")
            exit()


def write_default_config():
    default_config = {
        "db_name": "FACET_db",
        "outdir_name": "FACET_out"
    }

    try:
        with open(parent_dir+config_filename, 'w') as config:
            config.write(json.dumps(default_config))
            config.write("\n")

    except FileExistsError:
        print("%s%s already exists! Delete the current file if you want to rewrite it." % (parent_dir, config_filename))


def try_config_key(config_dict, key):
    try:
        config_dict[key]
        return True
    except KeyError:
        return False
