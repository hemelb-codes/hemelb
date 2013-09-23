//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#include "lb/iolets/InOutLetVelocityAware.h"
#include "configuration/SimConfig.h"
#include "net/NetworkTopology.h"
#include "geometry/LatticeData.h"
#include "lb/MacroscopicPropertyCache.h"

namespace hemelb
{
  namespace lb
  {
    namespace iolets
    {

      void InOutLetVelocityAware::InitialiseNeighbouringSites(
                                                              geometry::neighbouring::NeighbouringDataManager *manager,
                                                              geometry::LatticeData* latDat,
                                                              hemelb::lb::MacroscopicPropertyCache* globalPropertyCache,
                                                              std::vector<site_t> invertedBoundaryList)
      {
        NDM = manager;
        latticeData = latDat;

        hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("InitialiseNeighbouringSites: invBList size = %i",
                                                                              invertedBoundaryList.size());

        sitesWhichNeighbourThisBoundary = invertedBoundaryList;

        /* We retrieve the local velocities from the PropertyCache. */
        propertyCache = globalPropertyCache;

        // For each site neighbouring this boundary...
        for (std::vector<site_t>::iterator site_iterator = sitesWhichNeighbourThisBoundary.begin(); site_iterator
            != sitesWhichNeighbourThisBoundary.end(); site_iterator++)
        {
          // ... find its home proc.
          hemelb::util::Vector3D < site_t > gc;

          latDat->GetGlobalCoordsFromGlobalNoncontiguousSiteId(*site_iterator, gc);
          const proc_t neighbourSiteHomeProc = latDat->GetProcIdFromGlobalCoords(gc);

          // If this is a real proc, register the site as needed with the neighbouring data manager.
          if (neighbourSiteHomeProc != BIG_NUMBER2 && neighbourSiteHomeProc
              != net::NetworkTopology::Instance()->GetLocalRank())
          {
            manager->RegisterNeededSite(*site_iterator,
                                        geometry::neighbouring::RequiredSiteInformation(true));
          }
        }
      }

      InOutLetVelocityAware::InOutLetVelocityAware() :
        InOutLetMultiscale(), NDM(NULL), propertyCache(NULL), sitesWhichNeighbourThisBoundary(),
            latticeData(NULL)

      {
      }

      InOutLetVelocityAware::InOutLetVelocityAware(const InOutLetVelocityAware &other) :
        InOutLetMultiscale(other), NDM(other.NDM), propertyCache(other.propertyCache),
            latticeData(other.latticeData)
      {
      }

      InOutLetVelocityAware::~InOutLetVelocityAware()
      {
      }

      InOutLet* InOutLetVelocityAware::Clone() const
      {
        InOutLetVelocityAware* copy = new InOutLetVelocityAware(*this);
        return copy;
      }

      PhysicalVelocity_deprecated InOutLetVelocityAware::InOutLetVelocityAware::GetVelocity() const
      {
        util::Vector3D<float> thisnormal;
        thisnormal = const_cast<InOutLetVelocityAware*> (this)->GetNormal();

        util::Vector3D < distribn_t > totalVelocity(0.);

        int MyProcNumber = net::NetworkTopology::Instance()->GetLocalRank();

        /* Apply CalcDensityAndVelocity to extract velocities and add them all up.
         * We're not (yet) using weights or normalisation here, or conversion to physical units */
        for (std::vector<site_t>::const_iterator site_iterator =
            sitesWhichNeighbourThisBoundary.begin(); site_iterator
            != sitesWhichNeighbourThisBoundary.end(); site_iterator++)
        {
          ///TODO: Get the velocities from sites on this core too... right now we're only getting neighbouring data.
          hemelb::util::Vector3D < site_t > gc;
          latticeData->GetGlobalCoordsFromGlobalNoncontiguousSiteId(*site_iterator, gc);

          const proc_t neighbourSiteHomeProc = latticeData->GetProcIdFromGlobalCoords(gc);

          if (neighbourSiteHomeProc != BIG_NUMBER2 && neighbourSiteHomeProc == MyProcNumber
              && MyProcNumber > 0)
          {
            const util::Vector3D<distribn_t>
                v =
                    propertyCache->velocityCache.Get(latticeData->GetLocalContiguousIdFromGlobalNoncontiguousId(*site_iterator));

            totalVelocity += v;
          }
          else if (neighbourSiteHomeProc != MyProcNumber)
          {
            distribn_t* f;
            f
                = latticeData->GetNeighbouringData().GetSite(*site_iterator).GetFOld(latticeData->GetLatticeInfo().GetNumVectors());

            for (Direction direction = 0; direction < latticeData->GetLatticeInfo().GetNumVectors(); ++direction)
            {
              // TODO This is currently getting momentum, not velocity. We'd
              // need to divide by density if we wanted to sum just the velocity.
              totalVelocity
                  += util::Vector3D<distribn_t>(latticeData->GetLatticeInfo().GetVector(direction))
                      * f[direction];
            }
          }
        }
        /* Dot product of the total velocity with the boundary normal. */
        // TODO Should this dot product not be taken inside the loop for a surface integral?
        return totalVelocity.Dot(normal);
      }
    }
  }
}
