#ifndef HEMELB_NET_MIXINS_STORINGNET_HPP
#define HEMELB_NET_MIXINS_STORINGNET_HPP
#include "net/BaseNet.h"
#include "log/Logger.h"
namespace hemelb
{
  namespace net
  {
    class StoringNet : public virtual BaseNet
    {
      public:

        template<class T>
        void RequestSend(T* pointer, int count, proc_t rank)
        {
          if (count > 0)
          {

            sendProcessorComms[rank].push_back(BaseRequest(pointer, count, MpiDataType<T>()));
          }
        }

        template<class T>
        void RequestReceive(T* pointer, int count, proc_t rank)
        {
          if (count > 0)
          {
            receiveProcessorComms[rank].push_back(BaseRequest(pointer, count, MpiDataType<T>()));
          }
        }

        template<class T>
        void RequestGatherVSend(T* buffer, int count, proc_t toRank)
        {
          gatherVSendProcessorComms[toRank].push_back(BaseRequest(buffer, count, MpiDataType<T>()));
        }

        template<class T>
        void RequestGatherReceive(T* buffer)
        {
          gatherReceiveProcessorComms.push_back(ScalarRequest(buffer, MpiDataType<T>()));
        }

        template<class T>
        void RequestGatherSend(T* buffer, proc_t toRank)
        {
          gatherSendProcessorComms[toRank].push_back(ScalarRequest(buffer, MpiDataType<T>()));
        }

        template<class T>
        void RequestGatherVReceive(T* buffer, int * displacements, int *counts)
        {
          gatherVReceiveProcessorComms.push_back(GatherVReceiveRequest(buffer,
                                                                       displacements,
                                                                       counts,
                                                                       MpiDataType<T>()));
        }
      protected:
        /**
         * Struct representing all that's needed to successfully communicate with another processor.
         */
        class BaseRequest
        {
          public:
            void * Pointer;
            int Count;
            MPI_Datatype Type;
            BaseRequest(void *pointer, int count, MPI_Datatype type) :
                Pointer(pointer), Count(count), Type(type)
            {
            }
        };

        class ScalarRequest
        {
          public:
            void * Pointer;
            MPI_Datatype Type;
            ScalarRequest(void *pointer, MPI_Datatype type) :
                Pointer(pointer), Type(type)
            {
            }
        };

        class GatherVReceiveRequest
        {
          public:
            void * Pointer;
            int * Displacements;
            int * Counts;
            MPI_Datatype Type;
            GatherVReceiveRequest(void *pointer, int *displacements, int *counts, MPI_Datatype type) :
                Pointer(pointer), Displacements(displacements), Counts(counts), Type(type)
            {
            }
        };

        template<class Request>
        class BaseProcComms : public std::vector<Request>
        {
          public:
            MPI_Datatype Type;
        };

        class ProcComms : public BaseProcComms<BaseRequest>
        {
          public:
            void CreateMPIType()
            {
              std::vector<MPI_Aint> displacements(size());
              std::vector<int> lengths;
              std::vector<MPI_Datatype> types;

              int lLocation = 0;

              MPI_Aint offset;
              MPI_Get_address(front().Pointer, &offset);

              for (iterator it = begin(); it != end(); ++it)
              {
                MPI_Get_address(it->Pointer, &displacements[lLocation]);
                displacements[lLocation] -= offset;
                ++lLocation;
                lengths.push_back(it->Count);
                types.push_back(it->Type);
              }
              // Create the type and commit it.
              MPI_Type_create_struct(this->size(), &lengths.front(), &displacements.front(), &types.front(), &Type);
              MPI_Type_commit(&Type);
            }
        };

        class GatherProcComms : public BaseProcComms<ScalarRequest>
        {

        };

        class GatherVReceiveProcComms : public BaseProcComms<GatherVReceiveRequest>
        {

        };

        std::map<proc_t, ProcComms> sendProcessorComms;
        std::map<proc_t, ProcComms> receiveProcessorComms;
        std::map<proc_t, ProcComms> gatherVSendProcessorComms;
        GatherVReceiveProcComms gatherVReceiveProcessorComms;
        std::map<proc_t, GatherProcComms> gatherSendProcessorComms;
        GatherProcComms gatherReceiveProcessorComms;

    };
  }
}
#endif
