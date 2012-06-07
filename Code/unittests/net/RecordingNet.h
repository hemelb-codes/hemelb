#ifndef HEMELB_UNITTESTS_NET_RECORDINGNET_H
#define HEMELB_UNITTESTS_NET_RECORDINGNET_H
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
            for (std::map<proc_t, ProcComms>::iterator it = receiveProcessorComms.begin();
                it != receiveProcessorComms.end(); ++it)
            {
              for (ProcComms::iterator message = it->second.begin(); message != it->second.end(); message++)
              {

                if (requiredReceipts[message->Rank].size() == 0)
                {
                  CPPUNIT_ASSERT(requiredReceipts[message->Rank].size() != 0);
                  return;
                }

                CPPUNIT_ASSERT(requiredReceipts[message->Rank].front().EnvelopeIdentical(*message));
                requiredReceipts[message->Rank].front().Unpack(*message);

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

          class LabelledRequest : public SimpleRequest
          {
            public:
              const std::string Label;
              LabelledRequest(void *pointer, int count, MPI_Datatype type, proc_t rank, const std::string &label) :
                SimpleRequest(pointer, count, type, rank), Label(label)
              {
              }
              virtual bool EnvelopeIdentical(const SimpleRequest & other)
              {
                bool this_ok = ( (Count == other.Count) && (Rank == other.Rank) && (Type == other.Type));
                if (!this_ok)
                {
                  std::cerr << "Envelope different: " << Label << " R: " << Rank << " C: " << Count << " T: "
                      << static_cast<void*>(Type) << " : " << " R" << other.Rank << " C " << other.Count << " T "
                      << static_cast<void*>(other.Type) << std::flush;
                }
                return this_ok;
              }
              virtual bool PayloadIdentical(const SimpleRequest & other)
              {
                // reduction
                bool ok = true;
                for (int element = 0; element < Count; element++)
                {
                  int size;
                  MPI_Type_size(other.Type, &size);
                  // The below use of unsigned char is not formally correct (due to the possibility of char not having alignment 1)
                  // But we cannot currently see a better solution to avoid compiler warnings from void* arithmetic.
                  bool this_ok = 0
                      == std::memcmp(static_cast<unsigned char*>(other.Pointer) + size * element,
                                     static_cast<unsigned char*>(Pointer) + size * element,
                                     size);
                  if (!this_ok)
                  {

                    std::cerr << "Unexpected data in request: " << Label << " R " << Rank << " C " << Count << " : "
                        << std::endl;
                    for (int i = 0; i < size; i++)
                    {
                      // The below use of unsigned char is not formally correct (due to the possibility of char not having alignment 1)
                      // But we cannot currently see a better solution to avoid compiler warnings from void* arithmetic.
                      std::cerr << i << " : "
                          << static_cast<int>(* (static_cast<unsigned char*>(Pointer) + size * element + i)) << " "
                          << static_cast<int>(* (static_cast<char*>(other.Pointer) + size * element + i)) << std::endl;
                    }
                    std::cerr << std::endl;
                  }
                  ok = ok && this_ok;
                }
                return ok;
              }
              virtual void Unpack(SimpleRequest & other)
              {
                for (int element = 0; element < Count; element++)
                {
                  int size;
                  MPI_Type_size(other.Type, &size);
                  // The below use of unsigned char is not formally correct (due to the possibility of char not having alignment 1)
                  // But we cannot currently see a better solution to avoid compiler warnings from void* arithmetic.
                  std::memcpy(static_cast<unsigned char*>(other.Pointer) + size * element,
                              static_cast<unsigned char*>(Pointer) + size * element,
                              size);
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
          NetMock(topology::Communicator & communicator) :
              BaseNet(communicator), RecordingNet()
          {
          }

      };
    }
  }
}
#endif
