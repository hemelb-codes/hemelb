#ifndef HEMELB_UNITTESTS_GEOMETRY_MOCKS_H
#define HEMELB_UNITTESTS_GEOMETRY_MOCKS_H
#include <vector>
#include <map>
#include <cstdlib>
#include <iostream>

#include "constants.h"
#include "mpiInclude.h"
#include "net/net.h"

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
            requiredReceipts[rank].push_back(BaseRequest(pointer, count, MpiDataType<T>(),rank));
          }
          template<class T> void RequireSend(T* pointer, unsigned int count, proc_t rank, const std::string &label = "")
          {
            requiredSends[rank].push_back(BaseRequest(pointer, count, MpiDataType<T>(),rank));
          }
          void ReceivePointToPoint()
          {

            for (std::map<proc_t, ProcComms>::iterator it = receiveProcessorComms.begin();
                it != receiveProcessorComms.end(); ++it)
            {
              for (ProcComms::iterator message = it->second.begin(); message != it->second.end(); message++)
              {
                CPPUNIT_ASSERT(requiredReceipts[message->Rank].size() != 0);
                if (requiredReceipts[message->Rank].size() == 0)
                {
                  return;
                }

                CPPUNIT_ASSERT(requiredReceipts[message->Rank].front().EnvelopeIdentical(*message));
                requiredReceipts[message->Rank].front().Unpack(*message);
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
                CPPUNIT_ASSERT(requiredSends[message->Rank].size() != 0);
                if (requiredSends[message->Rank].size() == 0)
                {
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
        private:
          std::map<proc_t, ProcComms> requiredReceipts;
          std::map<proc_t, ProcComms> requiredSends;
      }
      ;

      class NetMock : public InterfaceDelegationNet<RecordingNet>,
                      public GathersViaPointPoint
      {
        public:
          NetMock() :
             BaseNet(), RecordingNet()
          {
          }

      };
    }
  }
}
#endif
