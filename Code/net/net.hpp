#ifndef HEMELB_NET_NET_HPP
#define HEMELB_NET_NET_HPP
#include "net/net.h"
namespace hemelb
{
  namespace net
  {
    template<class T>
    void Net::RequestSend(T* pointer, int count, proc_t rank)
    {
      if (count > 0)
      {
        if (sendReceivePrepped)
        {
          std::cerr << "Error: tried to add send-data after the datatype was already constructed. This is a bug.\n";
          exit(1);
        }
        sendProcessorComms[rank].push_back(BaseRequest(pointer, count, MpiDataType<T>()));
      }
    }

    template<class T>
    void Net::RequestReceive(T* pointer, int count, proc_t rank)
    {
      if (count > 0)
      {
        if (sendReceivePrepped)
        {
          std::cerr << "Error: tried to add receive-data after the datatype was already constructed. This is a bug.\n";
          exit(1);
        }
        receiveProcessorComms[rank].push_back(BaseRequest(pointer, count, MpiDataType<T>()));
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
    void Net::RequestGatherVReceive(std::vector<std::vector<T> > &buffer)
    {
      std::vector<int> displacements;
      std::vector<int> counts;
      for (typename std::vector<std::vector<T> >::iterator buffer_iterator = buffer.begin();
          buffer_iterator != buffer.end(); buffer++)
      {
        displacements.push_back(&buffer_iterator->front() - &buffer.front());
        counts.push_back(buffer_iterator->size());
      }
      RequestGatherVReceive(&buffer.front(), &displacements.front(), &counts.front());
    }

    template<class T>
    void Net::RequestGatherReceive(std::vector<T> &buffer)
    {

      RequestGatherReceive(&buffer.front());
    }

    template<class T>
    void Net::RequestGatherSend(T& value, proc_t toRank)
    {
      RequestGatherSend(&value, toRank);
    }

    template<class T>
    void Net::RequestGatherVSend(std::vector<T> &payload, proc_t toRank)
    {
      RequestGatherVSend(&payload.front(), payload.size(), toRank);
    }

    template<class T>
    void Net::RequestGatherVSend(T* buffer, int count, proc_t toRank)
    {
      gatherVSendProcessorComms[toRank].push_back(BaseRequest(buffer, count, MpiDataType<T>()));
    }

    template<class T>
    void Net::RequestGatherReceive(T* buffer)
    {
      gatherReceiveProcessorComms.push_back(ScalarRequest(buffer, MpiDataType<T>()));
    }

    template<class T>
    void Net::RequestGatherSend(T* buffer, proc_t toRank)
    {
      gatherSendProcessorComms[toRank].push_back(ScalarRequest(buffer, MpiDataType<T>()));
    }

    template<class T>
    void Net::RequestGatherVReceive(T* buffer, int * displacements, int *counts)
    {
      gatherVReceiveProcessorComms.push_back(GatherVReceiveRequest(buffer, displacements, counts, MpiDataType<T>()));
    }

  }
}

#endif // HEMELB_NET_NET_HPP
