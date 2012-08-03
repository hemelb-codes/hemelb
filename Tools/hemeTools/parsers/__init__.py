# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 

"""Package for reading/writing HemeLB file formats.

Beginning August 2011, binary formats will be identified with magic numbers
http://en.wikipedia.org/wiki/File_format#Magic_number

Here we define the generic HemeLB binary file identifier. It should be
the first 4 bytes of every (binary) file used for IO. It will be then
followed by another 4 bytes identifying the particular type/version,
that number being terminated by the EOF character (0x04).
"""

import numpy as np

HemeLbMagicNumber = 0x686c6221
