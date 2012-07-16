// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_UNITTESTS_NET_PHASED_MOCKCONCERN_H
#define HEMELB_UNITTESTS_NET_PHASED_MOCKCONCERN_H
#include "net/phased/Concern.h"
namespace hemelb
{
  namespace unittests
  {
    namespace net
    {
      namespace phased
      {
        using namespace hemelb::net::phased;
        class MockConcern : public Concern
        {
          public:
            MockConcern(const std::string &name) :
                calls(), name(name)
            {
            }

            bool CallAction(int action)
            {
              // this is where a real concern would switch on the action, and call the appropriate method
              calls.push_back(action);
              return true;
            }

            std::vector<int> const & ActionsCalled() const
            {
              return calls;
            }

          private:
            std::vector<int> calls;
            std::string name;
        };
      }
    }
  }
}
#endif
