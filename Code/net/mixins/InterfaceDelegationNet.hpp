#ifndef HEMELB_NET_MIXINS_INTERFACEDELEGATIONNET_HPP
#define HEMELB_NET_MIXINS_INTERFACEDELEGATIONNET_HPP
namespace hemelb
{
  namespace net
  {
    /***
     * Define the external template-based interface to be used by client classes.
     * We define a series of templates with nice C++ style interfaces which delegate to C-style template interfaces
     * And then, delegate the templated C-style interfaces to nontemplated interfaces taking an MPI datatype.
     */
    template<class BaseNet> class InterfaceDelegationNet : public virtual BaseNet
    {
      public:
        template<class T>
        void RequestSendV(std::vector<T> &payload, proc_t toRank)
        {
          RequestSend(&payload[0], payload.size(), toRank);
        }

        template<class T>
        void RequestSendR(T& value, proc_t toRank)
        {
          RequestSend(&value, 1, toRank);
        }

        template<class T>
        void RequestReceiveR(T& value, proc_t fromRank)
        {
          RequestReceive(&value, 1, fromRank);
        }

        template<class T>
        void RequestReceiveV(std::vector<T> &payload, proc_t toRank)
        {
          RequestReceive(&payload[0], payload.size(), toRank);
        }

        template<class T>
        void RequestGatherVReceive(std::vector<std::vector<T> > &buffer)
        {
          std::vector<int> & displacements = this->GetDisplacementsBuffer();
          std::vector<int> &counts = this->GetCountsBuffer();

          for (typename std::vector<std::vector<T> >::iterator buffer_iterator = buffer.begin();
              buffer_iterator != buffer.end(); buffer_iterator++)
          {
            // Ensure each vector has some underlying array, even if it's unused.
            buffer_iterator->reserve(1);
            displacements.push_back(&buffer_iterator->front() - &buffer.front().front());

            counts.push_back(buffer_iterator->size());
          }
          RequestGatherVReceive(&buffer.front().front(), &displacements.front(), &counts.front());
        }

        template<class T>
        void RequestGatherReceive(std::vector<T> &buffer)
        {
          // Ensure vector has some underlying array, even if it's unused.
          buffer.reserve(1);
          RequestGatherReceive(&buffer.front());
        }

        template<class T>
        void RequestGatherSend(T& value, proc_t toRank)
        {
          RequestGatherSend(&value, toRank);
        }

        template<class T>
        void RequestGatherVSend(std::vector<T> &payload, proc_t toRank)
        {
          RequestGatherVSend(&payload.front(), payload.size(), toRank);
        }

        template<class T>
        void RequestSend(T* pointer, int count, proc_t rank)
        {
            BaseNet::RequestSend(pointer, count, rank, MpiDataType<T>());
        }

        template<class T>
        void RequestReceive(T* pointer, int count, proc_t rank)
        {
            BaseNet::RequestReceive(pointer, count, rank, MpiDataType<T>());
        }

        /*
         * Blocking gathers are implemented in MPI as a single call for both send/receive
         * But, here we separate send and receive parts, since this interface may one day be used for
         * nonblocking collectives.
         */

        template<class T>
        void RequestGatherVSend(T* buffer, int count, proc_t toRank)
        {
            BaseNet::RequestGatherVSend(buffer, count, toRank, MpiDataType<T>());
        }

        template<class T>
        void RequestGatherReceive(T* buffer)
        {
            BaseNet::RequestGatherReceive(buffer, MpiDataType<T>());
        }

        template<class T>
        void RequestGatherSend(T* buffer, proc_t toRank)
        {
            BaseNet::RequestGatherSend(buffer, toRank, MpiDataType<T>());
        }

        template<class T>
        void RequestGatherVReceive(T* buffer, int * displacements, int *counts)
        {
            BaseNet::RequestGatherVReceive(buffer, displacements, counts, MpiDataType<T>());
        }

    };
  }
}
#endif
