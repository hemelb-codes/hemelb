#!/usr/bin/env python
# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 

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
config=yaml.load(open(os.path.join(localroot,'Tools','analysis','config_defaults.yml')))
config.update(yaml.load(open(os.path.join(localroot,'Tools','analysis','config.yml'))))

dc=yaml.load(open(os.path.join(localroot,'Tools','analysis','logging.yml')))
dc['handlers']['parsing']['filename']=os.path.expanduser(os.path.join(config['reports_path'],'parsing.log'))