#ifndef HEMELB_NET_H
#define HEMELB_NET_H

#include "constants.h"
#include "mpiInclude.h"
#include "D3Q15.h"
#include "SimConfig.h"

#include "lb/GlobalLatticeData.h"
#include "lb/LocalLatticeData.h"
#include "topology/NetworkTopology.h"

class Net
{
  public:
    Net(hemelb::topology::NetworkTopology * iTopology);
    ~Net();

    void Initialise(hemelb::lb::GlobalLatticeData &iGlobLatDat,
                    hemelb::lb::LocalLatticeData* &bLocalLatDat);

    void
    ReceiveFromNeighbouringProcessors(hemelb::lb::LocalLatticeData &bLocalLatDat);
    void
    SendToNeighbouringProcessors(hemelb::lb::LocalLatticeData &bLocalLatDat);
    void
    UseDataFromNeighbouringProcs(hemelb::lb::LocalLatticeData &bLocalLatDat);

    int my_inner_sites;
    int my_inner_collisions[COLLISION_TYPES];
    int my_inter_collisions[COLLISION_TYPES];
    MPI_Status status[4];
    double bm_time;
  private:
    void GetThisRankSiteData(const hemelb::lb::GlobalLatticeData & iGlobLatDat,
                             unsigned int *& bThisRankSiteData);
    void InitialiseNeighbourLookup(hemelb::lb::LocalLatticeData *bLocalLatDat,
                                   short int **bSharedFLocationForEachProc,
                                   const unsigned int *iSiteDataForThisRank,
                                   const hemelb::lb::GlobalLatticeData & iGlobLatDat);
    void CountCollisionTypes(const hemelb::lb::GlobalLatticeData & iGlobLatDat,
                             const unsigned int *lThisRankSiteData);
    void InitialisePointToPointComms(short int **& lSharedFLocationForEachProc);

    int *f_recv_iv;
    int err;

    hemelb::SimConfig* mSimConfig;
    hemelb::topology::NetworkTopology * mNetworkTopology;

    MPI_Request **req;

};

#endif // HEMELB_NET_H
