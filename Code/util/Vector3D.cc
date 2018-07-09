
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <cstdlib>
#include <iostream>
#include "util/Vector3D.h"

namespace hemelb
{
  namespace util
  {

    void Vector3DBase::SetIndexErrorHandler(Vector3DBase::HandlerFunction* func)
    {
      Vector3DBase::handler = func;
    }

    void Vector3DBase::HandleIndexError(int direction)
    {
      if (Vector3DBase::handler == NULL)
      {
        std::cout << "Vector3D index error handler not set when index error "
          "occurred" << std::endl;
        std::exit(1);
      }
      else
      {
        Vector3DBase::handler(direction);
      }
    }

    namespace detail
    {
      bool CheckNextChar(std::istream& i, char c)
      {
        char got = 0;
        i >> got;
        if (got != c)
        {
          i.setstate(i.failbit);
          return false;
        }
        return true;
      }
    }
  }
}
