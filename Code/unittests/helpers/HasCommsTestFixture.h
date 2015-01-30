//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_HELPERS_HASCOMMSTESTFIXTURE_H
#define HEMELB_UNITTESTS_HELPERS_HASCOMMSTESTFIXTURE_H
#include <cppunit/TestFixture.h>
#include "net/IOCommunicator.h"

namespace hemelb
{
  namespace unittests
  {
    namespace helpers
    {
      class HasCommsTestFixture : public CppUnit::TestFixture
      {

        public:
          void setUp()
          {
            //hemelbCommunicator = net::IOCommunicator::Instance();
          }

          void tearDown()
          {
            //hemelbCommunicator = nullptr;
          }

          static void Init(const net::IOCommunicator& inst)
          {
            hemelbCommunicator = &inst;
          }

        protected:
          static const net::IOCommunicator& Comms()
          {
            return *hemelbCommunicator;
          }
        private:
          static const net::IOCommunicator* hemelbCommunicator;
      };

      const net::IOCommunicator* HasCommsTestFixture::hemelbCommunicator = nullptr;
    }
  }
}
#endif // HEMELB_UNITTESTS_HELPERS_HASCOMMSTESTFIXTURE_H
