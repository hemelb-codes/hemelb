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
#include "SimulationMaster.h"

namespace hemelb
{
  namespace multiscale
  {
    /***
     * Instead of adding multiscale functionality to the standard simulation master, we keep this here,
     * so the main code can be read without thinking about multiscale.
     */
    template<class Intercommunicator>
    class MultiscaleSimulationMaster : public SimulationMaster<>
    {
      public:
        MultiscaleSimulationMaster(hemelb::configuration::CommandLine &options,
                                   const net::IOCommunicator& ioComm,
                                   Intercommunicator & aintercomms) :
            SimulationMaster(options, ioComm), intercomms(aintercomms),
                multiscaleIoletType("inoutlet")
        {
          // We only have one shared object type so far, an iolet.
          lb::iolets::InOutLetMultiscale::DefineType(multiscaleIoletType);

          hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::OnePerCore>("CONSTRUCTOR: inlet and outlet count: %d and %d",
                                                                               inletValues->GetLocalIoletCount(),
                                                                               outletValues->GetLocalIoletCount());
          hemelb::log::Logger::Log<hemelb::log::Warning, hemelb::log::OnePerCore>("inlets: %d",
                                                                                  inletValues->GetLocalIolet(0)->IsCommsRequired(),
                                                                                  inletValues->GetLocalIolet(0)->GetDensityMax(),
                                                                                  inletValues->GetLocalIolet(0)->GetPressureMax());
          hemelb::log::Logger::Log<hemelb::log::Warning, hemelb::log::OnePerCore>("outlets: %d",
                                                                                  outletValues->GetLocalIolet(0)->IsCommsRequired(),
                                                                                  outletValues->GetLocalIolet(0)->GetDensityMax(),
                                                                                  outletValues->GetLocalIolet(0)->GetPressureMax());

          // we only want to register those iolets which are needed on this process.
          // Fortunately, the BoundaryValues instance has worked this out for us.
          for (unsigned int i = 0; i < inletValues->GetLocalIoletCount(); i++)
          {
            // could be a if dynamic_cast<> rather than using a castable? virtual method pattern, if we prefer.
            if (inletValues->GetLocalIolet(i)->IsRegistrationRequired())
            {
              static_cast<lb::iolets::InOutLetMultiscale*>(inletValues->GetLocalIolet(i))->Register(intercomms,
                                                                                                    multiscaleIoletType);
            }
          }
          for (unsigned int i = 0; i < outletValues->GetLocalIoletCount(); i++)
          {
            if (outletValues->GetLocalIolet(i)->IsRegistrationRequired())
            {
              static_cast<lb::iolets::InOutLetMultiscale*>(outletValues->GetLocalIolet(i))->Register(intercomms,
                                                                                                     multiscaleIoletType);
            }
          }
          // We only have one shared object type so far, an iolet.
          lb::iolets::InOutLetMultiscale::DefineType(multiscaleIoletType);

          /* Process 0 has a list of all the Iolets. The count of all this is highly useful to pre-size all the
           * needed arrays later on, so we are broadcasting this to all the other processes. */
          std::vector<unsigned> GlobalIoletCount;
          GlobalIoletCount.push_back(inletValues->GetLocalIoletCount());
          GlobalIoletCount.push_back(outletValues->GetLocalIoletCount());
          ioComms.Broadcast(GlobalIoletCount, 0);

          std::vector<std::vector<site_t> > invertedInletBoundaryList(GlobalIoletCount[0]);
          std::vector<std::vector<site_t> > invertedOutletBoundaryList(GlobalIoletCount[1]);

          hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("inlets start %i/%i",
                                                                                inletValues->GetLocalIoletCount(),
                                                                                GlobalIoletCount[0]);
          hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("outlets start %i/%i",
                                                                                outletValues->GetLocalIoletCount(),
                                                                                GlobalIoletCount[1]);

          //TODO: Throw a warning when process 0 count mismatches with the aggregate of the others.
          bool velocity = false;
          if (velocity == true)
          {
            /* Do not include the non-iolet adjacent sites (resp. MidFluid and Wall-adjacent). */
            long long int offset = latticeData->GetMidDomainCollisionCount(0)
                + latticeData->GetMidDomainCollisionCount(1);

            /* Do include iolet adjacent sites (inlet) */
            long long int ioletsSiteCount = latticeData->GetMidDomainCollisionCount(2);
            invertedInletBoundaryList = PopulateInvertedBoundaryList(latticeData,
                                                                     invertedInletBoundaryList,
                                                                     offset,
                                                                     ioletsSiteCount);

            offset += latticeData->GetMidDomainCollisionCount(2);
            /* Do include iolet adjacent sites (outlet) */
            ioletsSiteCount = latticeData->GetMidDomainCollisionCount(3);
            invertedOutletBoundaryList = PopulateInvertedBoundaryList(latticeData,
                                                                      invertedOutletBoundaryList,
                                                                      offset,
                                                                      ioletsSiteCount);

            offset += latticeData->GetMidDomainCollisionCount(3);
            /* Do include iolet adjacent sites (inlet-wall) */
            ioletsSiteCount = latticeData->GetMidDomainCollisionCount(4);
            invertedInletBoundaryList = PopulateInvertedBoundaryList(latticeData,
                                                                     invertedInletBoundaryList,
                                                                     offset,
                                                                     ioletsSiteCount);

            offset += latticeData->GetMidDomainCollisionCount(4);
            /* Do include iolet adjacent sites (outlet-wall) */
            ioletsSiteCount = latticeData->GetMidDomainCollisionCount(5);
            invertedOutletBoundaryList = PopulateInvertedBoundaryList(latticeData,
                                                                      invertedOutletBoundaryList,
                                                                      offset,
                                                                      ioletsSiteCount);

            hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("Populated inlets (numinlets/sizeinlet0): %i/%i",
                                                                                  invertedInletBoundaryList.size(),
                                                                                  invertedInletBoundaryList[0].size());
            hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("Populated outlets (numoutlets/sizeoutlet0): %i/%i",
                                                                                  invertedOutletBoundaryList.size(),
                                                                                  invertedOutletBoundaryList[0].size());
          }

          //TODO: Debug
          //invertedInletBoundaryList =
          //    ExchangeAndCompleteInverseBoundaryList(invertedInletBoundaryList);
          //invertedOutletBoundaryList =
          //    ExchangeAndCompleteInverseBoundaryList(invertedOutletBoundaryList);

          //hemelb::lb::MacroscopicPropertyCache& propertyCache =
          //    latticeBoltzmannModel->GetPropertyCache();

          //WORKAROUND: Reserve space in the velocityCache to ensure that all elements are properly allocated.
          //If this is not set, then velocityCache will contain 0 elements at the start even though the
          //MacroscopicPropertyCache reports a number of site count higher than 0.
          //propertyCache->velocityCache.Reserve(propertyCache->GetSiteCount());

          // we only want to register those iolets which are needed on this process.
          // Fortunately, the BoundaryValues instance has worked this out for us.
          for (unsigned int i = 0; i < inletValues->GetLocalIoletCount(); i++)
          {
            hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("1) %i %i",
                                                                                  i,
                                                                                  GlobalIoletCount[0]);
            // could be a if dynamic_cast<> rather than using a castable? virtual method pattern, if we prefer.
            if (inletValues->GetLocalIolet(i)->IsRegistrationRequired())
            {
              hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("2) inlets: %i %i %i",
                                                                                    invertedInletBoundaryList.size(),
                                                                                    invertedInletBoundaryList[0].size(),
                                                                                    i);
              static_cast<lb::iolets::InOutLetMultiscale*>(inletValues->GetLocalIolet(i))->Register(intercomms,
                                                                                                    multiscaleIoletType);
              hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("3) inlets: %i %i",
                                                                                    invertedInletBoundaryList.size(),
                                                                                    invertedInletBoundaryList[i].size());
              /*static_cast<lb::iolets::InOutLetVelocityAware*>(inletValues->GetLocalIolet(i))->InitialiseNeighbouringSites(neighbouringDataManager,
               latticeData,
               &propertyCache,
               invertedInletBoundaryList[i]);*/
            }
          }

          for (unsigned int i = 0; i < outletValues->GetLocalIoletCount(); i++)
          {
            hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("1) %i %i",
                                                                                  i,
                                                                                  GlobalIoletCount[1]);
            if (outletValues->GetLocalIolet(i)->IsRegistrationRequired())
            {
              hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("2) outlets: %i %i %i",
                                                                                    invertedOutletBoundaryList.size(),
                                                                                    invertedOutletBoundaryList[0].size(),
                                                                                    i);
              static_cast<lb::iolets::InOutLetMultiscale*>(outletValues->GetLocalIolet(i))->Register(intercomms,
                                                                                                     multiscaleIoletType);
              hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("3) outlets: %i %i",
                                                                                    invertedOutletBoundaryList.size(),
                                                                                    invertedOutletBoundaryList[i].size());
              /*static_cast<lb::iolets::InOutLetVelocityAware*>(outletValues->GetLocalIolet(i))->InitialiseNeighbouringSites(neighbouringDataManager,
               latticeData,
               &propertyCache,
               invertedOutletBoundaryList[i]);*/
            }
          }

          hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("MSMaster ShareICs started...");
          intercomms.ShareInitialConditions();
          hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("MSMaster Init finished!");
        }

        void PrintVectorList(std::vector<std::vector<site_t> > v)
        {
          hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::OnePerCore>("Printing Vector List:");
          for (unsigned int i = 0; i < v.size(); i++)
          {
            for (unsigned int j = 0; j < v[i].size(); j++)
            {
              hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::OnePerCore>("Boundary: %i % i %i",
                                                                                   i,
                                                                                   j,
                                                                                   v[i][j]);
            }
          }
        }

        void DoTimeStep()
        {
          bool advance = intercomms.DoMultiscale(GetState()->GetTime());
          hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("At time step %i, should advance %i, time %f",
                                                                              GetState()->GetTimeStep(),
                                                                              static_cast<int>(advance),
                                                                              GetState()->GetTime());

          if (advance)
          {
            /* NOTE: Following triggers DoComms in each InOutLetMultiscale.
             * IoLetMS is aggressive with this, and DoComms will actually
             * complete all the communications, not just initiate them.
             * This is to prevent any inconsistent state in the coupling
             * (it's hard enough to get the physics right with a consistent
             * state ;)). */

            hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("inlet and outlet count: %d and %d",
                                                                                  inletValues->GetLocalIoletCount(),
                                                                                  outletValues->GetLocalIoletCount());
            hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("inlets: %d",
                                                                                  inletValues->GetLocalIolet(0)->IsCommsRequired(),
                                                                                  inletValues->GetLocalIolet(0)->GetDensityMax(),
                                                                                  inletValues->GetLocalIolet(0)->GetPressureMax());
            hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("outlets: %d",
                                                                                  outletValues->GetLocalIolet(0)->IsCommsRequired(),
                                                                                  outletValues->GetLocalIolet(0)->GetDensityMax(),
                                                                                  outletValues->GetLocalIolet(0)->GetPressureMax());

            SetCommsRequired(inletValues, true);
            SetCommsRequired(outletValues, true);

            inletValues->RequestComms();
            outletValues->RequestComms();
            SetCommsRequired(inletValues, false);
            SetCommsRequired(outletValues, false);

            for (unsigned int i = 0; i < inletValues->GetLocalIoletCount(); i++)
            {
              hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("Inlet[%i]: Measured Density is %f. Pressure is %f.",
                                                                                    i,
                                                                                    inletValues->GetLocalIolet(i)->GetDensity(GetState()->GetTimeStep()),
                                                                                    inletValues->GetLocalIolet(i)->GetPressureMax());
            }
            for (unsigned int i = 0; i < outletValues->GetLocalIoletCount(); i++)
            {
              hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("Outlet[%i]: Measured Density is %f. Pressure is %f.",
                                                                                    i,
                                                                                    outletValues->GetLocalIolet(i)->GetDensity(GetState()->GetTimeStep()),
                                                                                    outletValues->GetLocalIolet(i)->GetPressureMax());
            }

            /* Temporary Orchestration hardcode for testing 1/100 step ratio
             * TODO: Make an orchestration system for the multiscale coupling. */
            //for (int i = 0; i < 100; i++)
            //{
            hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::Singleton>("Step: HemeLB advanced to time %f.",
                                                                                 GetState()->GetTime());
            SimulationMaster::DoTimeStep();
            //}
          }
          else
          {
            hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::Singleton>("HemeLB waiting pending multiscale siblings.");
            return;
          };
        }
        Intercommunicator &intercomms;
        typename Intercommunicator::IntercommunicandTypeT multiscaleIoletType;

      private:

        /* Loops over iolets to set the need for communications. */
        void SetCommsRequired(hemelb::lb::iolets::BoundaryValues* ioletValues, bool b)
        {
          hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("Starting SetCommsRequired.");
          for (unsigned int i = 0; i < ioletValues->GetLocalIoletCount(); i++)
          {
            hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("In loop: %d",
                                                                                  ioletValues->GetLocalIoletCount());
            hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("A: iolet %d %d",
                                                                                  i,
                                                                                  (ioletValues->GetLocalIolet(i))->IsCommsRequired());
            hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("B: iolet %d %d",
                                                                                  i,
                                                                                  static_cast<lb::iolets::InOutLetMultiscale*>(ioletValues->GetLocalIolet(i))->IsCommsRequired());
            dynamic_cast<lb::iolets::InOutLetMultiscale*>(ioletValues->GetLocalIolet(i))->SetCommsRequired(b);
            hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("done with SetCommsRequired iteration.");

          }
          hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("Finishing SetCommsRequired.");
        }

        //Populate an invertedBoundaryList
        //out: iBL
        //in: LatticeData, [Site Object]->SiteData->GetBoundaryID,
        //MPI_Gatherv if data is only this process.
        std::vector<std::vector<site_t> > PopulateInvertedBoundaryList(
            hemelb::geometry::LatticeData* latticeData,
            std::vector<std::vector<site_t> > invertedBoundaryList, int offset, int ioletsSiteCount)
        {
          for (int i = 0; i < ioletsSiteCount; i++)
          {
            /* 1. Obtain Boundary ID number. */
            hemelb::geometry::SiteData s = latticeData->GetSite(offset + i).GetSiteData();
            int boundaryID = s.GetIoletId();

            /* 2. Grow the list to an appropriate size if needed. */
            while ( ((int) invertedBoundaryList.size()) <= boundaryID)
            {
              hemelb::log::Logger::Log<hemelb::log::Warning, hemelb::log::OnePerCore>("WARNING: Growing the invertedBoundaryList, because we created in wrongly in MultiscaleSimulation Master.");
              std::vector<site_t> a(0);
              invertedBoundaryList.push_back(a);
              if (invertedBoundaryList.size() > 100000)
              {
                hemelb::log::Logger::Log<hemelb::log::Warning, hemelb::log::OnePerCore>("ERROR: invertedBoundaryList is growing to ridiculous proportions due to a faulty boundaryID.");
                exit(-1);
              }
            }

            /* 3. Insert this site in the inverted Boundary List. */
            invertedBoundaryList[boundaryID].push_back(latticeData->GetGlobalNoncontiguousSiteIdFromGlobalCoords(latticeData->GetSite(offset
                + i).GetGlobalSiteCoords()));
          }
          return invertedBoundaryList;
        }

        std::vector<std::vector<site_t> > ExchangeAndCompleteInverseBoundaryList(
            std::vector<std::vector<site_t> > inList)
        {
          std::vector<std::vector<site_t> > outList;
          int *recvSizes = new int[ioComms.Size()];
          int *recvDispls = new int[ioComms.Size()];

          /* TODO: ASSUMPTION:
           * inList.size() is equal everywhere. This is not necessarily the case.
           * Use an AllReduce MAX and resize inList accordingly to make the remnant
           * of the code work here for unequal inList sizes. */

          for (unsigned int i = 0; i < inList.size(); i++)
          {
            int sendSize = ((int) inList[i].size());
            site_t *sendList = new site_t[inList[i].size()];
            for (unsigned int j = 0; j < inList[i].size(); j++)
            {
              sendList[j] = inList[i][j];
            }
            HEMELB_MPI_CALL(MPI_Allgather,
                            ( &sendSize, 1, MPI_INT, recvSizes, 1, MPI_INT, ioComms ));

            int64_t totalSize = 0;

            int np = ioComms.Size();
            int64_t offset = 0;

            for (int j = 0; j < np; j++)
            {
              totalSize += recvSizes[j];
              recvDispls[j] = offset;
              offset += recvSizes[j];
            }

            site_t *recvList = new site_t[totalSize]; //inList[i].size()

            HEMELB_MPI_CALL(MPI_Allgatherv,
                            ( sendList, inList[i].size(), MPI_LONG_LONG, recvList, recvSizes, recvDispls, MPI_LONG_LONG, ioComms ));

            std::vector<site_t> subList;
            for (int j = 0; j < totalSize; j++)
            {
              subList.push_back(recvList[j]);
            }
            outList.push_back(subList);
          }

          return outList;
        }
    };

  }
}

#endif // HEMELB_MULTISCALE_MULTISCALE_SIMULATION_MASTER_H
