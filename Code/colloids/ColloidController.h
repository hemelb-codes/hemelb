#ifndef HEMELB_COLLOIDS_COLLOIDCONTROLLER_H
#define HEMELB_COLLOIDS_COLLOIDCONTROLLER_H

#include <vector>
#include "net/net.h"
#include "geometry/LatticeData.h"
#include "util/Vector3D.h"
#include "geometry/ReadResult.h"

namespace hemelb
{
  namespace colloids
  {
    class ColloidController // : public net::IteratedAction
    {
      public:
        /** constructor - currently only initialises the neighbour list */
        ColloidController(net::Net* net,
                          geometry::LatticeData* latDatLBM,
                          geometry::GeometryReadResult* gmyResult);

        /** destructor - releases resources allocated by this class */
        ~ColloidController();

      private:
        /** enables simplified general point-to-point communication via MPI */
        net::Net* net;

        /** holds fluid information for local sites, i.e. the velocity distribution values */
        geometry::LatticeData* latDat;

        /** cached copy of local rank (obtained from topology) */
        proc_t localRank;

        /** maximum separation from a colloid of sites used in its fluid velocity interpolation */
        const static site_t REGION_OF_INFLUENCE = (site_t)2;

        /** a vector of the processors that might be interested in
            particles near the edge of this processor's sub-domain */
        std::vector<proc_t> neighbourProcessors;

        /** a list of relative 3D vectors that defines the sites within a region of influence */
        typedef std::vector<util::Vector3D<site_t> > Neighbourhood;

        /** obtains the neighbourhood for a particular region of influence defined by distance */
        Neighbourhood GetNeighbourhoodVectors(site_t distance);

        /** determines the list of neighbour processors
            i.e. processors that are within the region of influence of the local domain's edge
            i.e. processors that own at least one site in the neighbourhood of a local site */
        void InitialiseNeighbourList(geometry::GeometryReadResult* gmyResult,
                                     Neighbourhood neighbourhood);

        /** get local coordinates and the owner rank for a site from its global coordinates */
        bool GetLocalInformationForGlobalSite(geometry::GeometryReadResult* gmyResult,
                                              util::Vector3D<site_t> globalLocationForSite,
                                              site_t* blockIdForSite,
                                              site_t* localSiteIdForSite,
                                              proc_t* ownerRankForSite);
    };
  }
}

#endif /* HEMELB_COLLOIDS_COLLOIDCONTROLLER */
