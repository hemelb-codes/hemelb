# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

add_library(hemelb_multiscale
  mpwide/MPWideIntercommunicator.cc
  )
target_link_libraries(hemelb_multiscale PUBLIC MPWide::MPWide)
