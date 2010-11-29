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

  int lProcCount = 128;
  int lMachineCount = 1;

  lNetTop.NeighbourIndexFromProcRank = new short int[1];
  lNetTop.ProcCountOnEachMachine = new int[1];
  lNetTop.MachineIdOfEachProc = new int[lProcCount];
  lNetTop.ProcessorCount = lProcCount;
  lNetTop.MachineCount = 1;

  MakeAGlobLatDat(lGlobLatDat);

  int lFluidSites = lGlobLatDat.BlockCount
      * lGlobLatDat.SitesPerBlockVolumeUnit;

  double lExpectedPerRank = ((double) lFluidSites) / ((double) lProcCount - 1);

  lTopMan.DecomposeDomain(lFluidSites, lNetTop, lGlobLatDat);

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
  oGlobLatDat.SetBasicDetails(8, 8, 8, 4);

  for (int x = 0; x < oGlobLatDat.GetXBlockCount(); x++)
  {
    for (int y = 0; y < oGlobLatDat.GetYBlockCount(); y++)
    {
      for (int z = 0; z < oGlobLatDat.GetZBlockCount(); z++)
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
