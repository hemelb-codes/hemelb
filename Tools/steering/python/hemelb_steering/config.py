# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

"""
Parse the configs defined in config.yml and config_user.yml
"""

import os
import yaml

path=os.path.dirname(os.path.abspath(__file__))
config=yaml.load(open(os.path.join(path, 'config.yml')))
config_user=yaml.load(open(os.path.join(path, 'config_user.yml')))
