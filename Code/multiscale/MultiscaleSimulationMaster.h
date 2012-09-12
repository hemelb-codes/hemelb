// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_MULTISCALE_MULTISCALESIMULATIONMASTER_H
#define HEMELB_MULTISCALE_MULTISCALESIMULATIONMASTER_H
#include <vector>
#include "multiscale/Intercommunicator.h"
#include "lb/boundaries/iolets/InOutLetVelocityAware.h"
#include "SimulationMaster.h"

/*Temporary addition */
#include <mpi.h>

namespace hemelb
{
  namespace multiscale
  {
    /***
     * Instead of adding multiscale functionality to the standard simulation master, we keep this here,
     * so the main code can be read without thinking about multiscale.
     */
    template<class Intercommunicator> class MultiscaleSimulationMaster : public SimulationMaster
    {
      public:
        MultiscaleSimulationMaster(hemelb::configuration::CommandLine &options, Intercommunicator & aintercomms) :
            SimulationMaster(options), intercomms(aintercomms), multiscaleIoletType("inoutlet")
        {
          // We only have one shared object type so far, an iolet.

          lb::boundaries::iolets::InOutLetVelocityAware::DefineType(multiscaleIoletType);

          std::vector<std::vector<site_t> > invertedInletBoundaryList;
          std::vector<std::vector<site_t> > invertedOutletBoundaryList;

          /* Do not include the non-iolet adjacent sites (resp. MidFluid and Wall-adjacent). */
          long long int offset = latticeData->GetMidDomainCollisionCount(0) + latticeData->GetMidDomainCollisionCount(1);
          /* Do include iolet adjacent sites (inlet) */
          long long int ioletsSiteCount = latticeData->GetMidDomainCollisionCount(2);
          invertedInletBoundaryList = PopulateInvertedBoundaryList(latticeData, invertedInletBoundaryList, offset, ioletsSiteCount);

          offset += latticeData->GetMidDomainCollisionCount(2);
          /* Do include iolet adjacent sites (outlet) */
          ioletsSiteCount = latticeData->GetMidDomainCollisionCount(3);
          invertedOutletBoundaryList = PopulateInvertedBoundaryList(latticeData, invertedOutletBoundaryList, offset, ioletsSiteCount);

          offset += latticeData->GetMidDomainCollisionCount(3);
          /* Do include iolet adjacent sites (inlet-wall) */
          ioletsSiteCount = latticeData->GetMidDomainCollisionCount(4);
          invertedInletBoundaryList = PopulateInvertedBoundaryList(latticeData, invertedInletBoundaryList, offset, ioletsSiteCount);

          offset += latticeData->GetMidDomainCollisionCount(4);
          /* Do include iolet adjacent sites (outlet-wall) */
          ioletsSiteCount = latticeData->GetMidDomainCollisionCount(5);
          invertedOutletBoundaryList = PopulateInvertedBoundaryList(latticeData, invertedOutletBoundaryList, offset, ioletsSiteCount);

          invertedInletBoundaryList  = ExchangeAndCompleteInverseBoundaryList(invertedInletBoundaryList);
          invertedOutletBoundaryList = ExchangeAndCompleteInverseBoundaryList(invertedOutletBoundaryList);

          // we only want to register those iolets which are needed on this process.
          // Fortunately, the BoundaryValues instance has worked this out for us.
          for (unsigned int i = 0; i < inletValues->GetLocalIoletCount(); i++)
          {
            // could be a if dynamic_cast<> rather than using a castable? virtual method pattern, if we prefer.
            if (inletValues->GetLocalIolet(i)->IsRegistrationRequired())
            {
              static_cast<lb::boundaries::iolets::InOutLetMultiscale*>(inletValues->GetLocalIolet(i))->Register(intercomms,
                                                                                                                multiscaleIoletType);

              static_cast<lb::boundaries::iolets::InOutLetVelocityAware*>(inletValues->GetLocalIolet(i))->InitialiseNeighbouringSites(neighbouringDataManager,
                                                                                                                                      latticeData,
                                                                                                                                      invertedInletBoundaryList[i]);
            }
          }

          for (unsigned int i = 0; i < outletValues->GetLocalIoletCount(); i++)
          {
            if (outletValues->GetLocalIolet(i)->IsRegistrationRequired())
            {
              static_cast<lb::boundaries::iolets::InOutLetMultiscale*>(outletValues->GetLocalIolet(i))->Register(intercomms,
                                                                                                                 multiscaleIoletType);
              static_cast<lb::boundaries::iolets::InOutLetVelocityAware*>(outletValues->GetLocalIolet(i))->InitialiseNeighbouringSites(neighbouringDataManager,
                                                                                                                                       latticeData,
                                                                                                                                       invertedOutletBoundaryList[i]);
            }
          }

          intercomms.ShareInitialConditions();
        }

        void DoTimeStep()
        {

          if (intercomms.DoMultiscale(GetState()->GetTime()))
          {
            SimulationMaster::DoTimeStep(); //This one hangs!
            hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("HemeLB advanced to time %f.",
                                                                                GetState()->GetTime());
          }
          else
          {
            hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("HemeLB waiting pending multiscale siblings.");
            return;
          };
        }
        Intercommunicator &intercomms;
        typename Intercommunicator::IntercommunicandTypeT multiscaleIoletType;

      private:
        std::vector<std::vector<site_t> > PopulateInvertedBoundaryList(hemelb::geometry::LatticeData* latticeData,
                                                                       std::vector<std::vector<site_t> > invertedBoundaryList,
                                                                       int offset, int ioletsSiteCount)
        {
          //Populate an invertedBoundaryList
          //out: iBL
          //in: LatticeData, [Site Object]->SiteData->GetBoundaryID,
          //MPI_Gatherv if data is only this process.

          for(int i=0; i<ioletsSiteCount; i++) {
            /* 1. Obtain Boundary ID number. */
            hemelb::geometry::SiteData s = latticeData->GetSite(offset+i).GetSiteData();
            int boundaryID = s.GetBoundaryId();

            /* 2. Grow the list to an appropriate size if needed. */
            while(invertedBoundaryList.size() <= boundaryID)
            {
              std::vector<site_t> a(0);
              invertedBoundaryList.push_back(a);
              if(invertedBoundaryList.size()>100000)
              {
                std::cerr << "ERROR: invertedBoundaryList is growing to ridiculous proportions due to a faulty boundaryID." << std::endl;
                exit(-1);
              }
            }

            /* 3. Insert this site in the inverted Boundary List. */
            invertedBoundaryList[boundaryID].push_back(i+offset);
          }
          return invertedBoundaryList;
        }

        std::vector<std::vector<site_t> > ExchangeAndCompleteInverseBoundaryList(std::vector<std::vector<site_t> >  inList) {

          std::vector<std::vector<site_t> > outList;

          int sendSize = inList.size();
          int *recvSizes = new int[hemelb::topology::NetworkTopology::Instance()->GetProcessorCount()];
          int *recvDispls = new int[hemelb::topology::NetworkTopology::Instance()->GetProcessorCount()];

          for(int i=0; i<inList.size(); i++) {

            site_t *sendList = new site_t[inList[i].size()];
            for(int j=0; j<inList[i].size(); j++) {
              sendList[j] = inList[i][j];
            }

            MPI_Allgather(&sendSize, 1, MPI_INT, recvSizes, 1, MPI_INT, hemelb::topology::NetworkTopology::Instance()->GetComms().GetCommunicator());

            int64_t totalSize = 0;

            int np = 0;
            MPI_Comm_size(hemelb::topology::NetworkTopology::Instance()->GetComms().GetCommunicator(), &np);
            int64_t offset = 0;

            for(int j=0; j<np; j++) {
              totalSize += recvSizes[i];
              recvDispls[i] = offset;
              offset += recvSizes[i];
            }

            site_t *recvList = new site_t[inList[i].size()];

            MPI_Allgatherv(sendList, inList.size(), MPI_LONG_LONG, recvList, recvSizes, recvDispls,
                           MPI_LONG_LONG, hemelb::topology::NetworkTopology::Instance()->GetComms().GetCommunicator());

            std::vector<site_t> subList(totalSize);
            for(int j=0; j<totalSize; j++) {
              subList[j] = recvList[j]; //I could have used push_back here.
            }

            outList.push_back(subList);
          }
          return outList;
        }
    };
  }
}

#endif // HEMELB_MULTISCALE_MULTISCALE_SIMULATION_MASTER_H
