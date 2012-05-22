#ifndef HEMELB_NET_NET_HPP
#define HEMELB_NET_NET_HPP
#include "net/net.h"
namespace hemelb
{
  namespace net
  {
    template<class T>
    void Net::RequestSend(T* oPointer, int iCount, proc_t iToRank)
    {
      if (iCount > 0)
      {
        if (sendReceivePrepped)
        {
          std::cerr << "Error: tried to add send-data after the datatype was already constructed. This is a bug.\n";
          exit(1);
        }

        ProcComms *lComms = GetProcComms(iToRank, true);

        AddToList(oPointer, iCount, lComms);
      }
    }

    template<class T>
    void Net::RequestReceive(T* oPointer, int iCount, proc_t iFromRank)
    {
      if (iCount > 0)
      {
        if (sendReceivePrepped)
        {
          std::cerr << "Error: tried to add receive-data after the datatype was already constructed. This is a bug.\n";
          exit(1);
        }

        ProcComms *lComms = GetProcComms(iFromRank, false);

        AddToList(oPointer, iCount, lComms);
      }
    }

    template<class T>
    void Net::RequestSend(std::vector<T> &payload, proc_t toRank)
    {
      RequestSend(&payload[0], payload.size(), toRank);
    }

    template<class T>
    void Net::RequestSend(T& value, proc_t toRank)
    {
      RequestSend(&value, 1, toRank);
    }

    template<class T>
    void Net::RequestReceive(T& value, proc_t fromRank)
    {
      RequestReceive(&value, 1, fromRank);
    }

    template<class T>
    void Net::RequestReceive(std::vector<T> &payload, proc_t toRank)
    {
      RequestReceive(&payload[0], payload.size(), toRank);
    }

    template<typename T>
    void Net::AddToList(T* dataToAdd, int dataLength, ProcComms *procCommsObjectToAddTo)
    {
      procCommsObjectToAddTo->PointerList.push_back(dataToAdd);
      procCommsObjectToAddTo->LengthList.push_back(dataLength);
      procCommsObjectToAddTo->TypeList.push_back(MpiDataType<T>());
    }
  }
}

#endif // HEMELB_NET_NET_HPP
