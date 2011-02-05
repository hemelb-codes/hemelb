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

        void Initialise(hemelb::lb::GlobalLatticeData &iGlobLatDat,
                        hemelb::lb::LocalLatticeData* &bLocalLatDat);

        void
        ReceiveFromNeighbouringProcessors(hemelb::lb::LocalLatticeData &bLocalLatDat);
        void
        SendToNeighbouringProcessors(hemelb::lb::LocalLatticeData &bLocalLatDat);
        void
        UseDataFromNeighbouringProcs(hemelb::lb::LocalLatticeData &bLocalLatDat);

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
          ProcessorComms *lComms = GetProcComms(iToRank);

          AddToList(oPointer, iCount, lComms->mSendPointerList, lComms->mSendLengthList,
                    lComms->mSendTypeList);
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
          ProcessorComms *lComms = GetProcComms(iFromRank);

          AddToList(oPointer, iCount, lComms->mReceivePointerList, lComms->mReceiveLengthList,
                    lComms->mReceiveTypeList);
        }

        void InitialiseSendReceive(hemelb::lb::LocalLatticeData &bLocalLatDat);

      private:
        /**
         * Struct representing all that's needed to successfully communicate with another processor.
         */
        struct ProcessorComms
        {
          public:
            std::vector<void*> mSendPointerList;
            std::vector<int> mSendLengthList;
            std::vector<MPI_Datatype> mSendTypeList;
            std::vector<void*> mReceivePointerList;
            std::vector<int> mReceiveLengthList;
            std::vector<MPI_Datatype> mReceiveTypeList;

            MPI_Datatype SendType;
            MPI_Datatype ReceiveType;
        };

        ProcessorComms* GetProcComms(int iRank);

        void AddToList(int* iNew,
                       int iLength,
                       std::vector<void*> &iPointerList,
                       std::vector<int> &iLengthList,
                       std::vector<MPI_Datatype> &iTypeList);

        void AddToList(double* iNew,
                       int iLength,
                       std::vector<void*> &iPointerList,
                       std::vector<int> &iLengthList,
                       std::vector<MPI_Datatype> &iTypeList);

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

        void EnsurePreparedToSendReceive();

        void CreateMPIType(std::vector<void*> &iPointerList,
                           std::vector<int> &iLengthList,
                           std::vector<MPI_Datatype> &iTypeList,
                           MPI_Datatype &oNewDatatype);

        bool preppedToSendReceive;

        int *f_recv_iv;
        int err;

        std::map<int, ProcessorComms*> mProcessorComms;

        hemelb::topology::NetworkTopology * mNetworkTopology;

        MPI_Request **req;
        MPI_Status status[4];
    };

  }
}

#endif // HEMELB_NET_H
