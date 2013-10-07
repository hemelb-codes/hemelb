// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

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
