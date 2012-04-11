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
      private:
        // enables simplified general point-to-point communication via MPI
        net::Net* mNet;

        // holds fluid information for local sites, i.e. the velocity distribution values
        geometry::LatticeData* mLatDat;

        proc_t mLocalRank;

        // a vector of the processors that might be interested in
        // particles near the edge of this processor's sub-domain
        std::vector<proc_t> mNeighbourProcessors;

        typedef std::vector<util::Vector3D<site_t> > neighbourhood_t;
        neighbourhood_t GetNeighbourhoodVectors(site_t distance);

        void InitialiseNeighbourList(geometry::GeometryReadResult* gmyResult,
                                     neighbourhood_t neighbourhood);

        bool GetLocalInformationForGlobalSite(geometry::GeometryReadResult* gmyResult,
                                              util::Vector3D<site_t> globalLocationForSite,
                                              site_t* blockIdForSite,
                                              site_t* localSiteIdForSite,
                                              proc_t* ownerRankForSite);

      public:
        // constructor - called by SimulationMaster::Initialise()
        ColloidController(net::Net* net,
                          geometry::LatticeData* latDatLBM,
                          geometry::GeometryReadResult* gmyResult);

        // destructor - called by SimulationMaster::~SimulationMaster()
        ~ColloidController();
    };
  }
}

#endif /* HEMELB_COLLOIDS_COLLOIDCONTROLLER */
