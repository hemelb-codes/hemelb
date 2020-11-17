# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

configure_file("Logger.in.cc" "${CMAKE_CURRENT_BINARY_DIR}/Logger.cc")
add_library(hemelb_log ${CMAKE_CURRENT_BINARY_DIR}/Logger.cc ${logger})
