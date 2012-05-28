"""
Parse the configs defined in config.yml and config_user.yml
"""

import os
import yaml
path=os.path.dirname(os.path.abspath(__file__))
config=yaml.load(open(os.path.join(path,'config.yml')))
config_user=yaml.load(open(os.path.join(path,'config_user.yml')))
