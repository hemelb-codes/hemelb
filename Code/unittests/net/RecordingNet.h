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
            requiredReceipts[rank].push_back(LabelledRequest(pointer, count, MpiDataType<T>(), rank, label));
          }
          template<class T> void RequireSend(T* pointer, unsigned int count, proc_t rank, const std::string &label = "")
          {
            requiredSends[rank].push_back(LabelledRequest(pointer, count, MpiDataType<T>(), rank, label));
          }
          void ReceivePointToPoint()
          {

            std::cerr << "Unloading receipts from buffer" << std::endl;
            for (std::map<proc_t, ProcComms>::iterator it = receiveProcessorComms.begin();
                it != receiveProcessorComms.end(); ++it)
            {
              std::cerr << "Unloading receipts from " << it->first << " from buffer" << std::endl;
              for (ProcComms::iterator message = it->second.begin(); message != it->second.end(); message++)
              {

                if (requiredReceipts[message->Rank].size() == 0)
                {
                  std::cerr << "Request received to receive from " << message->Rank << "But no messages expected"
                      << std::endl;
                  CPPUNIT_ASSERT(requiredReceipts[message->Rank].size() != 0);
                  return;
                }

                CPPUNIT_ASSERT(requiredReceipts[message->Rank].front().EnvelopeIdentical(*message));
                requiredReceipts[message->Rank].front().Unpack(*message);
                std::cerr << "Unloading receipt from buffer" << std::endl;
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
                  std::cerr << "Request received to send to " << message->Rank << "But no messages expected"
                      << std::endl;
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
        private:

          class LabelledRequest : public BaseRequest
          {
            public:
              const std::string & Label;
              LabelledRequest(void *pointer, int count, MPI_Datatype type, proc_t rank, const std::string &label) :
                  BaseRequest(pointer, count, type, rank), Label(label)
              {
              }
              virtual bool EnvelopeIdentical(const BaseRequest & other)
              {
                std::cerr << "Check envelope" << std::endl;
                bool this_ok = ( (Count == other.Count) && (Rank == other.Rank) && (Type == other.Type));
                if (!this_ok)
                {
                  std::cerr << "Envelope different: " << " R: " << Rank  << " C: " << Count << " T: " << static_cast<void*>(Type) << " : "
                      << " R" << other.Rank << " C " << other.Count << " T " << static_cast<void*>(other.Type) << std::flush;
                }
                return this_ok;
              }
              virtual bool PayloadIdentical(const BaseRequest & other)
              {
                // reduction
                bool ok = true;
                for (unsigned int element = 0; element < Count; element++)
                {
                  int size;
                  MPI_Type_size(other.Type, &size);
                  bool this_ok = 0==std::memcmp(other.Pointer + size * element, Pointer + size * element, size);
                  if (!this_ok)
                  {
                    std::cerr << "Unexpected data in request: " << " R " << Rank << " C " << Count << " : "
                        << std::endl;
                    for (unsigned int i = 0; i < size; i++)
                    {
                      std::cerr << i << " : " << static_cast<int>(*static_cast<char*>(Pointer + size * element + i)) << " "
                          << static_cast<int>(*static_cast<char*>(other.Pointer + size * element + i)) << std::endl;
                    }
                    std::cerr << std::endl;
                  }
                  ok = ok && this_ok;
                }
                return ok;
              }
              virtual void Unpack(BaseRequest & other)
              {
                for (unsigned int element = 0; element < Count; element++)
                {
                  int size;
                  MPI_Type_size(other.Type, &size);
                  std::memcpy(other.Pointer + size * element, Pointer + size * element, size);
                }
              }
          };
          std::map<proc_t, BaseProcComms<LabelledRequest> > requiredReceipts;
          std::map<proc_t, BaseProcComms<LabelledRequest> > requiredSends;
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
