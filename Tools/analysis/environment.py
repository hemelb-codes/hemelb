#!/usr/bin/env python
# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 
import yaml
import os
import logging
import logging.config


localroot = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
#Load and invoke the default non-machine specific config JSON dictionaries.
defaults_file = open(os.path.join(localroot, 'Tools', 'analysis', 'config_defaults.yml'))
config = yaml.load(defaults_file)
defaults_file.close()

user_file = open(os.path.join(localroot, 'Tools', 'analysis', 'config.yml'))
config.update(yaml.load(user_file))
user_file.close()

logging_settings_file = open(os.path.join(localroot, 'Tools', 'analysis', 'logging.yml'))
dc = yaml.load(logging_settings_file)
logging_settings_file.close()
dc['handlers']['parsing']['filename'] = os.path.expanduser(os.path.join(config['reports_path'], 'parsing.log'))
logging.config.dictConfig(dc)
