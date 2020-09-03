# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

# TODO: when cmake min version increases to 3.13, can drop ${CMAKE_CURRENT_SOURCE_DIR}
target_sources(hemelb-tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/main.cc)
add_subdirectory(helpers)
target_sources(hemelb-tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/SimulationMasterTests.cc)
add_subdirectory(configuration)
add_subdirectory(extraction)
add_subdirectory(io)
