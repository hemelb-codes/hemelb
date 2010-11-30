#ifndef HEMELB_NET_H
#define HEMELB_NET_H

#include "constants.h"
#include "mpiInclude.h"
#include "D3Q15.h"
#include "SimConfig.h"

#include "lb/GlobalLatticeData.h"
#include "lb/LocalLatticeData.h"
#include "topology/TopologyManager.h"
#include "topology/NetworkTopology.h"

class Net
{
  public:
    Net(hemelb::topology::NetworkTopology * iTopology,
        int &iArgumentCount,
        char* iArgumentList[]);
    ~Net();

    void Abort();

    void Initialise(hemelb::topology::NetworkTopology &iNetTop,
                    hemelb::lb::GlobalLatticeData &iGlobLatDat,
                    hemelb::lb::LocalLatticeData* &bLocalLatDat);

    void
        ReceiveFromNeighbouringProcessors(hemelb::lb::LocalLatticeData &bLocalLatDat);
    void
    SendToNeighbouringProcessors(hemelb::lb::LocalLatticeData &bLocalLatDat);
    void
    UseDataFromNeighbouringProcs(hemelb::lb::LocalLatticeData &bLocalLatDat);

    int err;
    int my_inner_sites; // Site on this process that do and do not need
    // information from neighbouring processors.
    int my_inner_collisions[COLLISION_TYPES]; // Number of collisions that only use data on this rank.
    int my_inter_collisions[COLLISION_TYPES]; // Number of collisions that require information from
    // other processors.

    unsigned int GetCollisionType(unsigned int site_data);
    MPI_Status status[4];
    double bm_time, fr_time;

  private:
    // NeighProc is part of the Net (defined later in this file).  This object is an element of an array
    // (called neigh_proc[]) and comprises information about the neighbouring processes to this process.

    int *f_recv_iv;
    int my_inter_sites;

    hemelb::SimConfig* mSimConfig;
    hemelb::topology::NetworkTopology * mNetworkTopology;

    MPI_Request **req;

};

#endif // HEMELB_NET_H
