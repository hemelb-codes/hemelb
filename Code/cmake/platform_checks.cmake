# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
include_guard()

include(gnu_bug)
include(mountain_lion_scandir)
include(intel_cpp11)

include(CheckCXXSourceCompiles)
include(CheckCXXSymbolExists)

CHECK_CXX_SOURCE_COMPILES("#include <cmath>\n int main(int c,char** v){ return isnan(1.0); }" HAVE_ISNAN)
CHECK_CXX_SOURCE_COMPILES("#include <cmath>\n int main(int c,char** v){ return std::isnan(1.0); }" HAVE_STD_ISNAN)
CHECK_CXX_SOURCE_COMPILES("#include <sys/time.h>\n#include <sys/resource.h>\nint main(int c,char** v){ rusage usage;\ngetrusage(RUSAGE_SELF, &usage);\nreturn usage.ru_maxrss; }" HAVE_RUSAGE)

CHECK_CXX_SOURCE_COMPILES("
#include <stdint.h>
#include <rpc/types.h>
#include <rpc/xdr.h>
int main(int count, char** v){
  char buffer[15] = \"aaaaaaaaaaaaa\";
  XDR xdr;
  xdrmem_create(&xdr, buffer, 32, XDR_ENCODE);
  uint16_t a;
  uint32_t b;
  uint64_t c;
  xdr_uint16_t(&xdr, &a);
  xdr_uint32_t(&xdr, &b);
  xdr_uint64_t(&xdr, &c);
  return b;
}" HAVE_XDRUINTXX_T)

