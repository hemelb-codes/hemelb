#include "tests.h"

#include <math.h>
#include <cstdlib>

#include "topology/TopologyManager.h"

int main(int argc, char *argv[])
{
  Tests::RunTests();
}

void Tests::RunTests()
{
  TestDomainDecomposition();
}

void Tests::TestDomainDecomposition()
{
  hemelb::topology::TopologyManager lTopMan;
  hemelb::topology::NetworkTopology lNetTop;
  hemelb::lb::GlobalLatticeData lGlobLatDat;
  hemelb::lb::LocalLatticeData lLocalLatDat;

  int lProcCount = 128;
  int lMachineCount = 1;

  lNetTop.NeighbourIndexFromProcRank = new short int[1];
  lNetTop.ProcCountOnEachMachine = new int[1];
  lNetTop.MachineIdOfEachProc = new int[lProcCount];
  lNetTop.ProcessorCount = lProcCount;
  lNetTop.MachineCount = 1;

  MakeAGlobLatDat(lGlobLatDat);

  int lFluidSites = lGlobLatDat.BlocksX * lGlobLatDat.BlocksY
      * lGlobLatDat.BlocksZ * lGlobLatDat.SitesPerBlockVolumeUnit;

  double lExpectedPerRank = ((double) lFluidSites) / ((double) lProcCount - 1);

  lTopMan.DecomposeDomain(lFluidSites, lNetTop, lGlobLatDat, lLocalLatDat);

  for (int ii = 1; ii < lProcCount; ii++)
  {
    if (std::abs(lExpectedPerRank
        - (double) lNetTop.FluidSitesOnEachProcessor[ii]) > 1.0)
    {
      printf(
             "Domain decomposition problem: Expecting %f sites per processor, but rank %d had %d.\n",
             lExpectedPerRank, ii, lNetTop.FluidSitesOnEachProcessor[ii]);
    }
  }
}

void Tests::MakeAGlobLatDat(hemelb::lb::GlobalLatticeData & oGlobLatDat)
{
  oGlobLatDat.BlocksX = oGlobLatDat.BlocksY = oGlobLatDat.BlocksZ = 8;
  oGlobLatDat.Blocks = new hemelb::lb::BlockData[8 * 8 * 8];
  oGlobLatDat.BlockSize = 4;
  oGlobLatDat.Log2BlockSize = 2;
  oGlobLatDat.SitesX = oGlobLatDat.SitesY = oGlobLatDat.SitesZ = 8 * 4;
  oGlobLatDat.SitesPerBlockVolumeUnit = oGlobLatDat.BlockSize
      * oGlobLatDat.BlockSize * oGlobLatDat.BlockSize;

  for (int x = 0; x < oGlobLatDat.BlocksX; x++)
  {
    for (int y = 0; y < oGlobLatDat.BlocksY; y++)
    {
      for (int z = 0; z < oGlobLatDat.BlocksZ; z++)
      {
        oGlobLatDat .Blocks[x * 64 + y * 8 + z].ProcessorRankForEachBlockSite
            = new int[oGlobLatDat.SitesPerBlockVolumeUnit];
        for (int lSite = 0; lSite < oGlobLatDat.SitesPerBlockVolumeUnit; lSite++)
        {
          oGlobLatDat .Blocks[x * 64 + y * 8 + z].ProcessorRankForEachBlockSite[lSite]
              = -1;
        }
      }
    }
  }
}
