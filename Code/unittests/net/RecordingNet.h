
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

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
      /**
       * A mock for the net class for testing.
       * You must first give a complete, per-rank-ordered list of the sends and
       * another of the receives that it is to carry out.
       *
       * For sends, it checks that the sent data is identical to that specified/
       *
       * For recvs, it delivers the data that you specified.
       *
       * For both, it checks metadata.
       *
       * Once you have completed a "round" of communications, call
       * ExpectationsAllCompleted to check that there are none outstanding.
       *
       */
      class RecordingNet : public virtual StoringNet
      {
        public:
          RecordingNet(const MpiCommunicator& comms) :
              BaseNet(comms), StoringNet(comms), requiredReceipts(), requiredSends()
          {
          }

          /**
           * Specify that this rank should receive a message
           * @param pointer - to the message data that will later be received
           * @param count - of number of elements in the message
           * @param rank - of the source of the message
           * @param label - used in error reporting (should be unique)
           */
          template<class T> void RequireReceive(T* pointer, unsigned int count, proc_t rank,
                                                const std::string &label = "")
          {
            requiredReceipts[rank].push_back(LabelledRequest(pointer,
                                                             count,
                                                             MpiDataType<T>(),
                                                             rank,
                                                             label));
          }
          /**
           * Specify that this rank should send a message
           * @param pointer - to the message data that should be sent
           * @param count - of the number of elements in the message
           * @param rank - of the destination of the message
           * @param label - used in error reporting (should be unique)
           */
          template<class T> void RequireSend(T* pointer, unsigned int count, proc_t rank,
                                             const std::string &label = "")
          {
            requiredSends[rank].push_back(LabelledRequest(pointer,
                                                          count,
                                                          MpiDataType<T>(),
                                                          rank,
                                                          label));
          }

          /**
           * Mock-execute queued receives.
           *
           * Does no actual communication, but copies in the mock data
           * supplied to RequireReceive. It also checks that the receives
           * match the required ones.
           */
          void ReceivePointToPoint()
          {
            for (std::map<proc_t, ProcComms>::iterator it = receiveProcessorComms.begin();
                it != receiveProcessorComms.end(); ++it)
            {
              for (ProcComms::iterator message = it->second.begin(); message != it->second.end();
                  message++)
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
          /**
           * Mock-execute queued sends.
           *
           * Does no actual communication, but checks that the sent data
           * matches the mock data supplied to RequireSend. It also checks that
           * the sends match the required ones.
           */
          void SendPointToPoint()
          {

            for (std::map<proc_t, ProcComms>::iterator it = sendProcessorComms.begin();
                it != sendProcessorComms.end(); ++it)
            {
              for (ProcComms::iterator message = it->second.begin(); message != it->second.end();
                  message++)
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

          /**
           * Mock-wait - clears the message queue
           */
          void WaitPointToPoint()
          {
            receiveProcessorComms.clear();
            sendProcessorComms.clear();
          }

          /**
           * Assert that all required sends and receives have occurred.
           */
          void ExpectationsAllCompleted()
          {
            for (std::map<proc_t, BaseProcComms<LabelledRequest> >::iterator receipts_from_core =
                requiredReceipts.begin(); receipts_from_core != requiredReceipts.end();
                receipts_from_core++)
            {
              CPPUNIT_ASSERT(0 == receipts_from_core->second.size());
            }
            for (std::map<proc_t, BaseProcComms<LabelledRequest> >::iterator sends_to_core =
                requiredSends.begin(); sends_to_core != requiredSends.end(); sends_to_core++)
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
