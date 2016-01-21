#!/usr/bin/env python
# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
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
