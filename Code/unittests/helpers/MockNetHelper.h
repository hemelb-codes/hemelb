// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_UNITTESTS_HELPERS_MOCKNETHELPER_H
#define HEMELB_UNITTESTS_HELPERS_MOCKNETHELPER_H

#include "unittests/net/NetMock.h"
namespace hemelb
{
  namespace unittests
  {
    namespace helpers
    {
      class MockNetHelper
      {
        protected:
          MockNetHelper():communicatorMock(NULL),netMock(NULL){}
          void setUp(const proc_t core_count, const proc_t current_core)
          {
            communicatorMock = new topology::Communicator(current_core, core_count);
            netMock = new net::NetMock(*communicatorMock);
          }
          void tearDown(){
            delete communicatorMock;
            delete netMock;
          }

          topology::Communicator *communicatorMock;
          net::NetMock *netMock;
      };

    }
  }
}

#endif // HEMELB_UNITTESTS_HELPERS_RANDOMSOURCE_H
