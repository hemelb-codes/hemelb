
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "net/mixins/gathers/SeparatedGathers.h"
namespace hemelb
{
  namespace net
  {
    SeparatedGathers::SeparatedGathers(const MpiCommunicator& comms) :
        BaseNet(comms), StoringNet(comms)
    {
    }

    void SeparatedGathers::WaitGathers()
    {
      /*
       * For each of the stored gather send requests...
       * Each of which, is all the requests to gather to a certain root...
       */
      for (std::map<proc_t, GatherProcComms>::iterator send_it = gatherSendProcessorComms.begin();
          send_it != gatherSendProcessorComms.end(); ++send_it)
      {
        /*
         * send_it->first is the rank to which the gather must be sent
         * send_it->second is a VECTOR of ScalarRequests, each a stored request defining Type and Pointer.
         * There is one entry, for each request to gather to that root.
         */
        /*
         * If I am sending to MYSELF, then I must assemble a gather request for sending and receiving myself
         */
        if (send_it->first == communicator.Rank())
        {
          /*
           * I may be the ROOT of several gathers.
           * If I am, then we hope, the sends and receives were defined in the same order.
           * So, we assume the first of the receive requests, matches the first of the send-to-self requests
           */
          int gather_index = 0;
          for (GatherProcComms::iterator receive_it = gatherReceiveProcessorComms.begin();
              receive_it != gatherReceiveProcessorComms.end(); ++receive_it)
          {
            ScalarRequest toself = send_it->second[gather_index];

            /*
             * So now, I can send/receive to myself, for this request
             */HEMELB_MPI_CALL(MPI_Gather,
                             (toself.Pointer, 1, toself.Type, receive_it->Pointer, 1, receive_it->Type, communicator.Rank(), communicator));
            ++gather_index;
          }
        }
        else
        /*
         * I am not sending to myself, so I just need to send each of the gathers to the given root
         * No attempt is made to coalesce several gathers to the same root.
         */
        {
          for (GatherProcComms::iterator req = send_it->second.begin();
              req != send_it->second.end(); req++)
          {
            HEMELB_MPI_CALL(MPI_Gather,
                            (req->Pointer, 1, req->Type, NULL, 1, req->Type, send_it->first, communicator));
          }
        }
      }

      gatherSendProcessorComms.clear();
      gatherReceiveProcessorComms.clear();
    }

    void SeparatedGathers::WaitGatherVs()
    {

      for (std::map<proc_t, ProcComms>::iterator send_it = gatherVSendProcessorComms.begin();
          send_it != gatherVSendProcessorComms.end(); ++send_it)
      {

        if (send_it->first == communicator.Rank())
        {
          int gather_index = 0;

          for (GatherVReceiveProcComms::iterator receive_it = gatherVReceiveProcessorComms.begin();
              receive_it != gatherVReceiveProcessorComms.end(); ++receive_it)
          {
            /***
             * Again, as for WaitGather, we assume that send/receive requests were placed in parallel ordering...
             */
            SimpleRequest toself = send_it->second[gather_index];

            HEMELB_MPI_CALL(MPI_Gatherv,
                            ( toself.Pointer, toself.Count, toself.Type, receive_it->Pointer, receive_it->Counts, receive_it->Displacements, receive_it->Type, communicator.Rank(), communicator));
            ++gather_index;
          }
        }
        else
        {

          for (ProcComms::iterator req = send_it->second.begin(); req != send_it->second.end();
              req++)
          {
            HEMELB_MPI_CALL(MPI_Gatherv,
                            ( req->Pointer, req->Count, req->Type, NULL, NULL, NULL, req->Type, send_it->first, communicator));
          }

        }
      }
      gatherVSendProcessorComms.clear();
      gatherVReceiveProcessorComms.clear();
    }
  }
}

