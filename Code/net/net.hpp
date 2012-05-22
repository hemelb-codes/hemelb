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

    template<class T>
    void Net::RequestGatherReceive(T* buffer, int * displacements, int *counts)
    {

    }

    template<class T>
    void Net::RequestGatherReceive(std::vector<std::vector<T> > &buffer)
    {
      std::vector<int> displacements;
      std::vector<int> counts;
      for (typename std::vector<std::vector<T> >::iterator buffer_iterator = buffer.begin(); buffer_iterator != buffer.end();
          buffer++)
      {
        displacements.push_back(&buffer_iterator->front() - &buffer.front());
        counts.push_back(buffer_iterator->size());
      }
      RequestGatherReceive(&buffer.front(), &displacements.front(), &counts.front());
    }

    template<class T>
    void Net::RequestGatherReceive(std::vector<T> &buffer)
    {
      std::vector<int> displacements;
      std::vector<int> counts;
      for (typename std::vector<std::vector<T> >::iterator buffer_iterator = buffer.begin(); buffer_iterator != buffer.end();
          buffer++)
      {
        displacements.push_back(&*buffer_iterator - &buffer->front());
        counts.push_back(1);
      }
      RequestGatherReceive(&buffer.front(), &displacements.front(), &counts.front());
    }

    template<class T>
    void Net::RequestGatherSend(T& value, proc_t toRank)
    {
      RequestGatherSend(&value,1,toRank);
    }

    template<class T>
    void Net::RequestGatherSend(std::vector<T> &payload, proc_t toRank)
    {
      RequestGatherSend(&payload.front(), payload.size(), toRank);
    }

    template<class T>
    void Net::RequestGatherSend(T* buffer, int count, proc_t toRank)
    {
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
