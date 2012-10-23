# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 

"""
Parse the configs defined in config.yml and config_user.yml
"""

import os
import yaml

path=os.path.dirname(os.path.abspath(__file__))
config=yaml.load(open(os.path.join(path, 'config.yml')))
config_user=yaml.load(open(os.path.join(path, 'config_user.yml')))
