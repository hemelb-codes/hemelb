#include "geometry/LatticeData.h"

namespace hemelb
{
  namespace geometry
  {
    LatticeData::LatticeData(const bool reserveSteeringCore,
                             int* totalFluidSites,
                             unsigned int siteMins[3],
                             unsigned int siteMaxes[3],
                             int* fluidSitePerProc,
                             lb::LbmParameters* bLbmParams,
                             SimConfig* bSimConfig,
                             double* lReadTime,
                             double* lDecomposeTime) :
      localLatDat(), globLatDat()

    {
      GeometryReader reader(reserveSteeringCore);

      reader.LoadAndDecompose(&globLatDat, totalFluidSites, siteMins, siteMaxes, fluidSitePerProc,
                              bLbmParams, bSimConfig, lReadTime, lDecomposeTime);

      int localRank;
      MPI_Comm_rank(MPI_COMM_WORLD, &localRank);

      localLatDat.Initialise(fluidSitePerProc[localRank]);
    }

    const double* LatticeData::GetNormalToWall(int iSiteIndex) const
    {
      return localLatDat.GetNormalToWall(iSiteIndex);
    }
  }
}
