#ifndef HEMELB_NET_H
#define HEMELB_NET_H

#include <vector>
#include <map>

#include "constants.h"
#include "mpiInclude.h"
#include "D3Q15.h"
#include "SimConfig.h"

#include "lb/GlobalLatticeData.h"
#include "lb/LocalLatticeData.h"
#include "topology/NetworkTopology.h"

namespace hemelb
{
  namespace net
  {

    class Net
    {
      public:
        Net(hemelb::topology::NetworkTopology * iTopology);
        ~Net();

        int* Initialise(hemelb::lb::GlobalLatticeData &iGlobLatDat,
                        hemelb::lb::LocalLatticeData* &bLocalLatDat);

        void Receive();
        void Send();
        void Wait(hemelb::lb::LocalLatticeData *bLocalLatDat);

        /**
         * Request that iCount entries of type T be included in the send to iToRank,
         * starting at oPointer.
         *
         * @param oPointer Pointer to first element to be sent.
         * @param iCount Number of element to be sent.
         * @param iToRank Rank to send to.
         */
        template<class T>
        void RequestSend(T* oPointer, int iCount, int iToRank)
        {
          if (sendReceivePrepped)
          {
            std::cerr
                << "Error: tried to add send-data after the datatype was already constructed. This is a bug.\n";
            exit(1);
          }

          ProcComms *lComms = GetProcComms(iToRank, true);

          AddToList(oPointer, iCount, lComms);
        }

        /**
         * Request that iCount entries of type T be included in the receive from iFromRank
         * starting at oPointer.
         *
         * @param oPointer Pointer to the first element of the array to receive into.
         * @param iCount Number of elements to receive into the array.
         * @param iFromRank Rank to receive from.
         */
        template<class T>
        void RequestReceive(T* oPointer, int iCount, int iFromRank)
        {
          if (sendReceivePrepped)
          {
            std::cerr
                << "Error: tried to add receive-data after the datatype was already constructed. This is a bug.\n";
            exit(1);
          }

          ProcComms *lComms = GetProcComms(iFromRank, false);

          AddToList(oPointer, iCount, lComms);
        }

      private:

        void GetThisRankSiteData(const hemelb::lb::GlobalLatticeData & iGlobLatDat,
                                 unsigned int *& bThisRankSiteData);
        void InitialiseNeighbourLookup(hemelb::lb::LocalLatticeData *bLocalLatDat,
                                       short int **bSharedFLocationForEachProc,
                                       const unsigned int *iSiteDataForThisRank,
                                       const hemelb::lb::GlobalLatticeData & iGlobLatDat);
        void CountCollisionTypes(hemelb::lb::LocalLatticeData * bLocalLatDat,
                                 const hemelb::lb::GlobalLatticeData & iGlobLatDat,
                                 const unsigned int * lThisRankSiteData);

        void InitialisePointToPointComms(short int **& lSharedFLocationForEachProc);

        /**
         * Struct representing all that's needed to successfully communicate with another processor.
         */
        struct ProcComms
        {
          public:
            std::vector<void*> PointerList;
            std::vector<int> LengthList;
            std::vector<MPI_Datatype> TypeList;

            void clear()
            {
              PointerList.clear();
              LengthList.clear();
              TypeList.clear();
            }

            MPI_Datatype Type;
        };

        ProcComms* GetProcComms(int iRank, bool iIsSend);

        void AddToList(int* iNew, int iLength, ProcComms* bMetaData);

        void AddToList(double* iNew, int iLength, ProcComms* bMetaData);

        void EnsurePreparedToSendReceive();

        void CreateMPIType(ProcComms *iMetaData);

        void EnsureEnoughRequests(unsigned int count);

        bool sendReceivePrepped;

        std::map<int, ProcComms*> mSendProcessorComms;
        std::map<int, ProcComms*> mReceiveProcessorComms;

        hemelb::topology::NetworkTopology * mNetworkTopology;

        // Requests and statuses available for general communication within the Net object (both
        // initialisation and during each iteration). Code using these must make sure
        // there are enough available. We do this in a way to minimise the number created
        // on each core, but also to minimise creation / deletion overheads.
        std::vector<MPI_Request> mRequests;
        std::vector<MPI_Status> mStatuses;
    };

  }
}

#endif // HEMELB_NET_H
