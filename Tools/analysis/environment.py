#!/usr/bin/env python
# encoding: utf-8
"""
environment.py

Created by James Hetherington on 2012-01-23.
Copyright (c) 2012 University College London. All rights reserved.
"""
import yaml
import os
import logging
import logging.config


localroot=os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
#Load and invoke the default non-machine specific config JSON dictionaries.
defaults_file = open(os.path.join(localroot,'Tools','analysis','config_defaults.yml'))
config=yaml.load(defaults_file)
defaults_file.close()

user_file = open(os.path.join(localroot,'Tools','analysis','config.yml'))
config.update(yaml.load(user_file))
user_file.close()

logging_settings_file = open(os.path.join(localroot,'Tools','analysis','logging.yml'))
dc=yaml.load(logging_settings_file)
logging_settings_file.close()
dc['handlers']['parsing']['filename']=os.path.expanduser(os.path.join(config['reports_path'],'parsing.log'))
logging.config.dictConfig(dc)
