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

namespace hemelb
{
  namespace lb
  {
    namespace boundaries
    {
      namespace iolets
      {

        PhysicalVelocity InOutLetVelocityAware::GetVelocity() const
        {
          util::Vector3D<float> thisnormal;
          thisnormal = const_cast<InOutLetVelocityAware*>(this)->GetNormal();

          distribn_t total_velocity[3];
          //          return 1.0;

          /* Apply CalcDensityAndVelocity to extract velocities and add them all up.
           * We're not (yet) using weights or normalisation here, or conversion to physical units */
          for (std::vector<site_t>::const_iterator site_iterator = sitesWhichNeighbourThisBoundary.begin();
              site_iterator != sitesWhichNeighbourThisBoundary.end(); site_iterator++)
          {

            distribn_t* f =
                latticeData->GetNeighbouringData().GetSite(*site_iterator).GetFOld(latticeData->GetLatticeInfo().GetNumVectors());
            //latticeData->CalculateDensityAndVelocity(f_data, &density, &velocity[0], &velocity[1], &velocity[2]);
            for (Direction direction = 0; direction < latticeData->GetLatticeInfo().GetNumVectors(); ++direction)
            {
              total_velocity[0] += latticeData->GetLatticeInfo().GetVector(direction)[0] * f[direction];
              total_velocity[1] += latticeData->GetLatticeInfo().GetVector(direction)[1] * f[direction];
              total_velocity[2] += latticeData->GetLatticeInfo().GetVector(direction)[2] * f[direction];
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
