
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_HELPERS_CPPUNITCOMPAREVECTORS_H
#define HEMELB_UNITTESTS_HELPERS_CPPUNITCOMPAREVECTORS_H
#include <cppunit/TestFixture.h>

// Let vectors be output to an ostream, so that CPPUNIT can assert equality on them

namespace CPPUNIT_NS
{
  template<class T> struct assertion_traits<std::vector<T> >
  {
      static bool equal(const std::vector<T>& x, const std::vector<T>& y)
      {
        return x == y;
      }

      // Note this vector print doesn't visually distinguish between ("1" "2" "3") and("1, 2" "3").
      static std::string toString(const std::vector<T>& values)
      {
        std::stringstream output;
        output << "[ " << std::flush;
        for (typename std::vector<T>::const_iterator value = values.begin(); value != values.end(); value++)
        {
          output << *value << ", " << std::flush;
        }
        output << "]" << std::flush;
        return output.str();
      }
  };
}
#endif
