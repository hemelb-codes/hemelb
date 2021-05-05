# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

"""Package for reading/writing HemeLB file formats.

Beginning August 2011, binary formats will be identified with magic numbers
http://en.wikipedia.org/wiki/File_format#Magic_number

Here we define the generic HemeLB binary file identifier. It should be
the first 4 bytes of every (binary) file used for IO. It will be then
followed by another 4 bytes identifying the particular type/version,
that number being terminated by the EOF character (0x04).
"""

import numpy as np

HemeLbMagicNumber = 0x686C6221
