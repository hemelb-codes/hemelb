
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_NET_MIXINS_INTERFACEDELEGATIONNET_H
#define HEMELB_NET_MIXINS_INTERFACEDELEGATIONNET_H
namespace hemelb
{
  namespace net
  {
    /***
     * Define the external template-based interface to be used by client classes.
     * We define a series of templates with nice C++ style interfaces which delegate to C-style template interfaces
     * And then, delegate the templated C-style interfaces to nontemplated interfaces taking an MPI datatype.
     */
    class InterfaceDelegationNet : public virtual BaseNet
    {
      public:
    InterfaceDelegationNet(comm::Communicator::ConstPtr comms) :
            BaseNet(comms)
        {
        }

        template<class T>
        void RequestSendV(const std::vector<T> &payload, proc_t toRank)
        {
          RequestSend(&payload[0], payload.size(), toRank);
        }

        template<class T>
        void RequestSendR(const T& value, proc_t toRank)
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
        void RequestSend(const T* pointer, int count, proc_t rank)
        {
          RequestSendImpl(const_cast<T*>(pointer), count, rank, comm::MpiDataType<T>());
        }

        template<class T>
        void RequestReceive(T* pointer, int count, proc_t rank)
        {
          RequestReceiveImpl(pointer, count, rank, comm::MpiDataType<T>());
        }

    };
  }
}
#endif
