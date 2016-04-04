
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_HELPERS_MOCKNETHELPER_H
#define HEMELB_UNITTESTS_HELPERS_MOCKNETHELPER_H

#include "unittests/net/NetMock.h"
namespace hemelb
{
  namespace unittests
  {
    namespace helpers
    {
      class MockMpiCommunicator : public net::MpiCommunicator
      {
        public:
          /***
           * Constructor for a dummy communicator
           * Can be useful for testing but can't actually be used
           * @param rank
           * @param size
           */
          MockMpiCommunicator(int rank_, int size_) :
            MpiCommunicator(), rank(rank_), size(size_)
          {

          }

          virtual inline int Rank() const
          {
            return rank;
          }
          virtual inline int Size() const
          {
            return size;
          }
        private:
          int rank, size;
      };

      class MockNetHelper
      {
        protected:
          MockNetHelper() :
            communicatorMock(NULL), netMock(NULL)
          {
          }
          void setUp(const proc_t core_count, const proc_t current_core)
          {
            communicatorMock = new MockMpiCommunicator(current_core, core_count);
            netMock = new net::NetMock(*communicatorMock);
          }
          void tearDown()
          {
            delete communicatorMock;
            delete netMock;
          }

          net::MpiCommunicator *communicatorMock;
          net::NetMock *netMock;
      };

    }
  }
}

#endif // HEMELB_UNITTESTS_HELPERS_RANDOMSOURCE_H
