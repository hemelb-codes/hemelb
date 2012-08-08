// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_COLLOIDS_COLLOIDCONTROLLER_H
#define HEMELB_COLLOIDS_COLLOIDCONTROLLER_H

#include <vector>
#include "net/net.h"
#include "net/IteratedAction.h"
#include "geometry/LatticeData.h"
#include "geometry/Geometry.h"
#include "io/xml/XmlAbstractionLayer.h"
#include "lb/MacroscopicPropertyCache.h"
#include "colloids/ParticleSet.h"
#include "util/Vector3D.h"

namespace hemelb
{
  namespace colloids
  {
    /** provides the control interface between colloid simulation and the rest of the system */
    class ColloidController : public net::IteratedAction
    {
      public:
        /** constructor - currently only initialises the neighbour list */
        ColloidController(const geometry::LatticeData& latDatLBM,
                          const lb::SimulationState& simulationState,
                          const geometry::Geometry& gmyResult,
                          io::xml::XmlAbstractionLayer& xml,
                          lb::MacroscopicPropertyCache& propertyCache);

        /** destructor - releases resources allocated by this class */
        ~ColloidController();

        /** overloaded from IteratedAction */
        void RequestComms();

        /** overloaded from IteratedAction */
        void EndIteration();

      private:
        /** cached copy of local rank (obtained from topology) */
        const proc_t localRank;

        /** holds the set of Particles that this processor knows about */
        ParticleSet* particleSet;

        /** maximum separation from a colloid of sites used in its fluid velocity interpolation */
        const static site_t REGION_OF_INFLUENCE = (site_t)2;

        /** a vector of the processors that might be interested in
            particles near the edge of this processor's sub-domain */
        std::vector<proc_t> neighbourProcessors;

        /** a list of relative 3D vectors that defines the sites within a region of influence */
        typedef std::vector<util::Vector3D<site_t> > Neighbourhood;

        /** obtains the neighbourhood for a particular region of influence defined by distance */
        const Neighbourhood GetNeighbourhoodVectors(site_t distance);

        /** determines the list of neighbour processors
            i.e. processors that are within the region of influence of the local domain's edge
            i.e. processors that own at least one site in the neighbourhood of a local site */
        void InitialiseNeighbourList(const geometry::LatticeData& latDatLBM,
                                     const geometry::Geometry& gmyResult,
                                     const Neighbourhood& neighbourhood);

        /** get local coordinates and the owner rank for a site from its global coordinates */
        bool GetLocalInformationForGlobalSite(const geometry::Geometry& gmyResult,
                                              const util::Vector3D<site_t>& globalLocationForSite,
                                              site_t* blockIdForSite,
                                              site_t* localSiteIdForSite,
                                              proc_t* ownerRankForSite);
    };
  }
}

#endif /* HEMELB_COLLOIDS_COLLOIDCONTROLLER */
