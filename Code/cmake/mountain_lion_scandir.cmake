# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
include_guard()

include(CheckCXXSourceCompiles)
# mountain lion changed the api of scandir from BSD to LINUX style
CHECK_CXX_SOURCE_COMPILES("
#include <string>
#include <dirent.h>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <unistd.h>
#include <sys/dir.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sstream>

typedef const struct direct direct_t;
int selectOnlyContents(direct_t *entry)
    {
    return 1;
    }
int main(int count, char** v)
    {
      std::string pathname = \"\\tmp\";
      struct direct **files;
      int file_count = scandir(pathname.c_str(), &files, selectOnlyContents, alphasort);
      return 0;
    }"
  LINUX_SCANDIR)
