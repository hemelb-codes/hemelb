#include "colloids/ColloidController.h"

namespace hemelb
{
  namespace colloids
  {

    // constructor - called by SimulationMaster::Initialise()
    ColloidController::ColloidController(net::Net* net,
                                         geometry::LatticeData* latDat,
                                         geometry::GeometryReadResult* gmyResult)
                                       : mNet(net), mLatDat(latDat)
    {
      // TODO: implementation - use the GeometryReadResult object to determine neighbour procs
    }

    // destructor
    ColloidController::~ColloidController()
    {
      /*
      mNet = NULL;
      mLatDat = NULL;
      */
    }
  };
}
