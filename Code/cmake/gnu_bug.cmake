# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
include_guard()

# Check for a serious compiler bug in GCC
# http://gcc.gnu.org/bugzilla/show_bug.cgi?id=50618

if(CMAKE_COMPILER_IS_GNUXX)
  include(CheckCXXSourceRuns)
  CHECK_CXX_SOURCE_RUNS("struct Base {
  const int text;
  Base() : text(1) {}
  Base(int aText) : text(aText) {}
};
struct SubA : public virtual Base {
protected:
  int x;
public:
  SubA(int aX) : x(aX) {}
};
class SubB : public virtual Base {};
struct Diamond : public SubA, public SubB {
    Diamond(int text) : Base(text), SubA(5), SubB() {}
    void printText() {
        if(text != 2)
          __builtin_abort();
        if(x!=5)
          __builtin_abort();
    }
};

int main(int, char**) {
    Diamond x(2);
    x.printText();
}" HAVE_GCC_WITHOUT_DIAMOND_BUG)
  
  if (NOT HAVE_GCC_WITHOUT_DIAMOND_BUG)
    message(SEND_ERROR "Your version of GCC has this bug http://gcc.gnu.org/bugzilla/show_bug.cgi?id=50618. HemeLB requires a working implementation of virtual base classes. Please use a different compiler.")
  endif()
endif()
