#ifndef HEMELB_COLLOIDS_COLLOIDCONTROLLER_H
#define HEMELB_COLLOIDS_COLLOIDCONTROLLER_H

#include "net/net.h"
#include "geometry/LatticeData.h"
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

      public:
        // constructor - called by SimulationMaster::Initialise()
        ColloidController(net::Net* net,
            geometry::LatticeData* latDat,
            geometry::GeometryReadResult* gmyResult);

        // destructor - called by SimulationMaster::~SimulationMaster()
        ~ColloidController();
    };

  }
}

#endif /* HEMELB_COLLOIDS_COLLOIDCONTROLLER */
