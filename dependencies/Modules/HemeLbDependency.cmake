# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
include_guard()

macro(hemelb_dependency NAME STEP)
  include(${HEMELB_DEPENDENCIES_PATH}/${NAME}/${STEP}.cmake)
endmacro()
