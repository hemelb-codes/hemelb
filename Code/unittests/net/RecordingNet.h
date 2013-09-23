// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_UNITTESTS_NET_RECORDINGNET_H
#define HEMELB_UNITTESTS_NET_RECORDINGNET_H
#include <vector>
#include <map>
#include <cstdlib>
#include <iostream>

#include "constants.h"
#include "net/mpi.h"
#include "net/net.h"
#include "unittests/net/LabelledRequest.h"

namespace hemelb
{
  namespace unittests
  {
    namespace net
    {
      using namespace hemelb::net;
      class RecordingNet : public virtual StoringNet
      {
        public:
          RecordingNet() :
              requiredReceipts(), requiredSends()
          {
          }

          template<class T> void RequireReceive(T* pointer, unsigned int count, proc_t rank, const std::string &label =
                                                    "")
          {
            requiredReceipts[rank].push_back(LabelledRequest(pointer, count, MpiDataType<T>(), rank, label));
          }
          template<class T> void RequireSend(T* pointer, unsigned int count, proc_t rank, const std::string &label = "")
          {
            requiredSends[rank].push_back(LabelledRequest(pointer, count, MpiDataType<T>(), rank, label));
          }
          void ReceivePointToPoint()
          {
            for (std::map<proc_t, ProcComms>::iterator it = receiveProcessorComms.begin();
                it != receiveProcessorComms.end(); ++it)
            {
              for (ProcComms::iterator message = it->second.begin(); message != it->second.end(); message++)
              {

                if (requiredReceipts[message->Rank].size() == 0)
                {
                  // It should not be the case, that a request is made, but there is no recorded assertion in the list of recorded assertions.
                  // i.e., we've popped off all the recorded requests, and there are none left to match the request we just received.
                  CPPUNIT_ASSERT(requiredReceipts[message->Rank].size() != 0);
                  return;
                }

                CPPUNIT_ASSERT(requiredReceipts[message->Rank].front().EnvelopeIdentical(*message));
                requiredReceipts[message->Rank].front().Unpack(*message);

                // This assertion just checks that "Unpack" has done it's job.
                // It shouldn't fail due to problems in tested code
                // Would only fail due to memory screwups or bad code in this class's Unpack method.
                CPPUNIT_ASSERT(requiredReceipts[message->Rank].front().PayloadIdentical(*message));

                requiredReceipts[message->Rank].pop_front();
              }

            }

          }
          void SendPointToPoint()
          {

            for (std::map<proc_t, ProcComms>::iterator it = sendProcessorComms.begin(); it != sendProcessorComms.end();
                ++it)
            {
              for (ProcComms::iterator message = it->second.begin(); message != it->second.end(); message++)
              {

                if (requiredSends[message->Rank].size() == 0)
                {
                  CPPUNIT_ASSERT(requiredSends[message->Rank].size() != 0);
                  return;
                }
                CPPUNIT_ASSERT(requiredSends[message->Rank].front().EnvelopeIdentical(*message));
                CPPUNIT_ASSERT(requiredSends[message->Rank].front().PayloadIdentical(*message));
                requiredSends[message->Rank].pop_front();
              }

            }
          }
          void WaitPointToPoint()
          {
            receiveProcessorComms.clear();
            sendProcessorComms.clear();
          }
          void ExpectationsAllCompleted()
          {
            for (std::map<proc_t, BaseProcComms<LabelledRequest> >::iterator receipts_from_core =
                requiredReceipts.begin(); receipts_from_core != requiredReceipts.end(); receipts_from_core++)
            {
              CPPUNIT_ASSERT(0 == receipts_from_core->second.size());
            }
            for (std::map<proc_t, BaseProcComms<LabelledRequest> >::iterator sends_to_core = requiredSends.begin();
                sends_to_core != requiredSends.end(); sends_to_core++)
            {
              CPPUNIT_ASSERT(0 == sends_to_core->second.size());
            }
          }
        private:


          std::map<proc_t, BaseProcComms<LabelledRequest> > requiredReceipts;
          std::map<proc_t, BaseProcComms<LabelledRequest> > requiredSends;
      }
      ;


    }
  }
}
#endif
