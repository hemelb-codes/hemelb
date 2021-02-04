# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
include_guard()

hemelb_dependency(ctemplate find)

if (CTEMPLATE_FOUND)
  message("Ctemplate already installed, no need to download")
  add_custom_target(ctemplate)
else()
  message("Ctemplate not installed, will build from source")
  find_file(CTEMPLATE_TARBALL ctemplate-2.3.tar.gz
    DOC "Path to download CTemplate (can be url http://)"
    PATHS ${HEMELB_DEPENDENCIES_PATH}/distributions
    )
  if(NOT CTEMPLATE_TARBALL)
    message("No ctemplate source found, will download.")
    set(CTEMPLATE_TARBALL http://github.com/OlafvdSpek/ctemplate/archive/ctemplate-2.3.tar.gz
      CACHE STRING "Path to download CTemplate (can be local file://)" FORCE)
  endif()
  
  set(CTEMPLATE_CONFIGURE_OPTIONS "" CACHE STRING "Extra configurations options for CTEMPLATE")
  option(CTEMPLATE_PATCH_VACOPY "Define va_copy macro through patch" OFF)
  option(CTEMPLATE_PATCH_ALIGN "Define GTL align macros as gnu" OFF)
  
  if (CTEMPLATE_PATCH_VACOPY)
    set(PATCH_COMMAND_VACOPY patch -p1 < ${HEMELB_DEPENDENCIES_PATH}/patches/ctemplate_vacopy.diff)
  else()
    set(PATCH_COMMAND_VACOPY echo novacopy)
  endif()
  if (CTEMPLATE_PATCH_ALIGN)
    set(PATCH_COMMAND_ALIGN patch -p1 < ${HEMELB_DEPENDENCIES_PATH}/patches/ctemplate_align.diff)
  else()
    set(PATCH_COMMAND_ALIGN echo noalign)
  endif()
  ExternalProject_Add(
    ctemplate
    INSTALL_DIR ${HEMELB_DEPENDENCIES_INSTALL_PATH}
    URL ${CTEMPLATE_TARBALL}
    CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR>  ${CTEMPLATE_CONFIGURE_OPTIONS}
    BUILD_COMMAND make -j${HEMELB_SUBPROJECT_MAKE_JOBS}
    BUILD_IN_SOURCE 1
    PATCH_COMMAND ${PATCH_COMMAND_ALIGN} && ${PATCH_COMMAND_VACOPY}
    )
endif()
