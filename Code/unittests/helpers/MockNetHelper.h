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
              communicatorMock(nullptr), netMock(nullptr)
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
