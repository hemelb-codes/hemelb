/*
 * InOutLetVelocityAware.cc
 *
 *  Created on: 3 Jul 2012
 *      Author: derek
 */

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

          std::cout << "invBList size: " << invBList.size() << std::endl;

          sitesWhichNeighbourThisBoundary = invBList;

          //std::cout << "sitesWhichNeighbourThisBoundary size: " << sitesWhichNeighbourThisBoundary.size() << std::endl;

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

          /* Apply CalcDensityAndVelocity to extract velocities and add them all up.
           * We're not (yet) using weights or normalisation here, or conversion to physical units */
          for (std::vector<site_t>::const_iterator site_iterator = sitesWhichNeighbourThisBoundary.begin();
              site_iterator != sitesWhichNeighbourThisBoundary.end(); site_iterator++)
          {

            ///TODO: Get the velocities from sites on this core too... right now we're only getting neighbouring data.
            hemelb::util::Vector3D<site_t> gc;
            latticeData->GetGlobalCoordsFromGlobalNoncontiguousSiteId(*site_iterator, gc);
            const proc_t neighbourSiteHomeProc = latticeData->GetProcIdFromGlobalCoords(gc);

            if (neighbourSiteHomeProc != BIG_NUMBER2
                && neighbourSiteHomeProc != topology::NetworkTopology::Instance()->GetLocalRank())
            {
              const util::Vector3D<distribn_t>& velocity = propertyCache->velocityCache.Get(*site_iterator);
              total_velocity[0] = velocity[0];
              total_velocity[1] = velocity[1];
              total_velocity[2] = velocity[2];
            }
            else
            {
              distribn_t* f =
                  latticeData->GetNeighbouringData().GetSite(*site_iterator).GetFOld(latticeData->GetLatticeInfo().GetNumVectors());
              //latticeData->CalculateDensityAndVelocity(f_data, &density, &velocity[0], &velocity[1], &velocity[2]);

              //std::cout << "IoletVA NumVectors: " << latticeData->GetLatticeInfo().GetNumVectors() << std::endl;

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
          }
          /* Dot product of the total velocity with the boundary normal. */
          // Should this dot product not be taken inside the loop for a surface integral?
          return (total_velocity[0] * normal[0]) + (total_velocity[1] * normal[1]) + (total_velocity[2] * normal[2]);
        }
      }
    }
  }
}
