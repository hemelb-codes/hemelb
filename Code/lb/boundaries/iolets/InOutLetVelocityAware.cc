//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#include "lb/boundaries/iolets/InOutLetVelocityAware.h"
#include "configuration/SimConfig.h"
#include "topology/NetworkTopology.h"
#include "geometry/LatticeData.h"
#include "lb/MacroscopicPropertyCache.h"

namespace hemelb
{
  namespace lb
  {
    namespace boundaries
    {
      namespace iolets
      {

        void InOutLetVelocityAware::InitialiseNeighbouringSites(geometry::neighbouring::NeighbouringDataManager *manager,
                                                                geometry::LatticeData* latDat,
                                                                hemelb::lb::MacroscopicPropertyCache* globalPropertyCache,
                                                                std::vector<site_t> invBList)
        {
          NDM = manager;
          latticeData = latDat;

          hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::OnePerCore>("InitialiseNeighbouringSites: invBList size = %i",
                                                                               invBList.size());

          sitesWhichNeighbourThisBoundary = invBList;

          /* We retrieve the local velocities from the PropertyCache. */
          propertyCache = globalPropertyCache;

          for (std::vector<site_t>::iterator site_iterator = sitesWhichNeighbourThisBoundary.begin();
              site_iterator != sitesWhichNeighbourThisBoundary.end(); site_iterator++)
          {
            hemelb::util::Vector3D<site_t> gc;

            latDat->GetGlobalCoordsFromGlobalNoncontiguousSiteId(*site_iterator, gc);
            const proc_t neighbourSiteHomeProc = latDat->GetProcIdFromGlobalCoords(gc);

            if (neighbourSiteHomeProc != BIG_NUMBER2
                && neighbourSiteHomeProc != topology::NetworkTopology::Instance()->GetLocalRank())
            {
              //hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::OnePerCore>("Site Needed (InitialiseNeighbouringSites)");
              manager->RegisterNeededSite(*site_iterator, geometry::neighbouring::RequiredSiteInformation(true));
            }
          }
        }

        PhysicalVelocity InOutLetVelocityAware::GetVelocity() const
        {
          util::Vector3D<float> thisnormal;
          thisnormal = const_cast<InOutLetVelocityAware*>(this)->GetNormal();

          distribn_t total_velocity[3] = { 0.0, 0.0, 0.0 };
          //          return 1.0;

          //std::cout << "IoletVA siteloop starts:" << std::endl;

          int MyProcNumber = topology::NetworkTopology::Instance()->GetLocalRank();

          int i = 0;

          /* Apply CalcDensityAndVelocity to extract velocities and add them all up.
           * We're not (yet) using weights or normalisation here, or conversion to physical units */
          for (std::vector<site_t>::const_iterator site_iterator = sitesWhichNeighbourThisBoundary.begin();
              site_iterator != sitesWhichNeighbourThisBoundary.end(); site_iterator++)
          {

            ///TODO: Get the velocities from sites on this core too... right now we're only getting neighbouring data.
            hemelb::util::Vector3D<site_t> gc;
            latticeData->GetGlobalCoordsFromGlobalNoncontiguousSiteId(*site_iterator, gc);

            const proc_t neighbourSiteHomeProc = latticeData->GetProcIdFromGlobalCoords(gc);

            //hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::OnePerCore>("Fluid site count: %i", latticeData->GetFluidSiteCountOnProc(1));

            if (neighbourSiteHomeProc != BIG_NUMBER2 && neighbourSiteHomeProc == MyProcNumber && MyProcNumber > 0)
            {
              /*hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::OnePerCore>("1. GetVel from PropertyCache for site %i, local id: %i, size: %i",
                                                                                   *site_iterator,
                                                                                   latticeData->GetLocalContiguousIdFromGlobalNoncontiguousId(*site_iterator),
                                                                                   propertyCache->GetSiteCount());
              */

              const util::Vector3D<distribn_t> v =
                  propertyCache->velocityCache.Get(latticeData->GetLocalContiguousIdFromGlobalNoncontiguousId(*site_iterator));

              /*hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::OnePerCore>("2. velocities are %f %f %f",
                                                                                   v[0],
                                                                                   v[1],
                                                                                   v[2]);*/
              total_velocity[0] += v[0];
              total_velocity[1] += v[1];
              total_velocity[2] += v[2];

              //hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::OnePerCore>("3. totalvel.%f", total_velocity[0]);
            }

            /*distribn_t* f;
             f = latticeData->GetSite(*site_iterator).GetFOld(latticeData->GetLatticeInfo().GetNumVectors()); */

            if (neighbourSiteHomeProc != MyProcNumber)
            {
              distribn_t* f;
              f =
                  latticeData->GetNeighbouringData().GetSite(*site_iterator).GetFOld(latticeData->GetLatticeInfo().GetNumVectors());

              //latticeData->CalculateDensityAndVelocity(f_data, &density, &velocity[0], &velocity[1], &velocity[2]);

              for (Direction direction = 0; direction < latticeData->GetLatticeInfo().GetNumVectors(); ++direction)
              {
                //std::cout << "Write out Vector directions and fs: " << std::endl;
                //std::cout << "Vectors: " << latticeData->GetLatticeInfo().GetVector(direction)[0] << " "
                //    << latticeData->GetLatticeInfo().GetVector(direction)[1] << " "
                //    << latticeData->GetLatticeInfo().GetVector(direction)[2] << " " << std::endl;
                //std::cout << "f[direction] = " << f[direction] << std::endl;
                total_velocity[0] += latticeData->GetLatticeInfo().GetVector(direction)[0] * f[direction];
                total_velocity[1] += latticeData->GetLatticeInfo().GetVector(direction)[1] * f[direction];
                total_velocity[2] += latticeData->GetLatticeInfo().GetVector(direction)[2] * f[direction];
              }
            }
            i++;
          }
          /* Dot product of the total velocity with the boundary normal. */
          // Should this dot product not be taken inside the loop for a surface integral?
          //hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::OnePerCore>("IOletVA GetVelocity(): Total velocity = %f", (total_velocity[0] * normal[0]) + (total_velocity[1] * normal[1]) + (total_velocity[2] * normal[2]));

          return (total_velocity[0] * normal[0]) + (total_velocity[1] * normal[1]) + (total_velocity[2] * normal[2]);
        }
      }
    }
  }
}
