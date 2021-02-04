# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
include_guard()

hemelb_dependency(cppunit find)
if (CPPUNIT_FOUND)
  message("CPPUnit already installed, no need to download")
  add_custom_target(cppunit)
else()
  message("CPPUnit not installed, will build from source")
  
  set(CPPUNIT_CONFIGURE_OPTIONS "" CACHE STRING "Extra configuration options for CPPUNIT")
  option(CPPUNIT_PATCH_LDL "Add -ldl option to cppunit" OFF)
  option(CPPUNIT_PATCH_DYNAMIC "Add -dynamic option to cppunit. Only works if CPPUNIT_PATCH_LDL is OFF." OFF)

  find_file(CPPUNIT_TARBALL cppunit-1.12.1.tar.gz 
    DOC "Path to download CPPUNIT (can be url http://)"
    PATHS ${HEMELB_DEPENDENCIES_PATH}/distributions
    )
  if(NOT CPPUNIT_TARBALL)
    message("No cppunit source found, will download.")
    set(CPPUNIT_TARBALL http://downloads.sourceforge.net/cppunit/cppunit-1.12.1.tar.gz
      CACHE STRING "Path to download CPPUNIT (can be local file://)" FORCE)
  endif()
  if(APPLE)
    set(RECONFIGURE_CPPUNIT_DEFAULT ON)
  else()
    set(RECONFIGURE_CPPUNIT_DEFAULT OFF)
  endif()
  option(RECONFIGURE_CPPUNIT "Include Reconfigure step for CPPUNIT" ${RECONFIGURE_CPPUNIT_DEFAULT})
  if (CPPUNIT_PATCH_LDL)
    set(PATCH_COMMAND_LDL patch -p1 < ${HEMELB_DEPENDENCIES_PATH}/patches/cppunit_ldl.diff)
  elseif (CPPUNIT_PATCH_DYNAMIC)
    set(PATCH_COMMAND_LDL patch -p1 < ${HEMELB_DEPENDENCIES_PATH}/patches/cppunit_dynamic.diff)
  else()
    set(PATCH_COMMAND_LDL echo noldl)
  endif()
  if (RECONFIGURE_CPPUNIT)
    set(PATCH_COMMAND_RECONFIGURE autoreconf -fvi) #autoreconf required on osx - based on contents of portfile)
  else()
    set(PATCH_COMMAND_RECONFIGURE echo noreconf)
  endif()
  ExternalProject_Add(
    cppunit
    INSTALL_DIR ${HEMELB_DEPENDENCIES_INSTALL_PATH}
    URL ${CPPUNIT_TARBALL}
    CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR> --mandir=<INSTALL_DIR>/share/man --disable-doxygen --disable-dot "${CPPUNIT_CONFIGURE_OPTIONS}"
    BUILD_IN_SOURCE 1
    BUILD_COMMAND make -j${HEMELB_SUBPROJECT_MAKE_JOBS}
    PATCH_COMMAND ${PATCH_COMMAND_LDL} && ${PATCH_COMMAND_RECONFIGURE}
    )
endif()
