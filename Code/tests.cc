#include "tests.h"

#include <cmath>
#include <cstdlib>

#include "topology/NetworkTopology.h"

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
  int lInputs = 1;
  char* lName = new char[6];
  strcpy(lName, "Tests");
  char ** lNamePtr = &lName;
  bool lSuccess;

  hemelb::topology::NetworkTopology * lNetTop = new hemelb::topology::NetworkTopology(&lInputs,
                                                                                      &lNamePtr,
                                                                                      &lSuccess);

  hemelb::lb::GlobalLatticeData lGlobLatDat;

  int lProcCount = lNetTop->GetProcessorCount();
  int lMachineCount = 1;

  lNetTop->NeighbourIndexFromProcRank = new short int[1];

  MakeAGlobLatDat(lGlobLatDat);

  int lFluidSites = lGlobLatDat.GetBlockCount() * lGlobLatDat.SitesPerBlockVolumeUnit;

  double lExpectedPerRank = ((double) lFluidSites) / ((double) lProcCount - 1);

  lNetTop->DecomposeDomain(lFluidSites, true, lGlobLatDat);

  for (int ii = 1; ii < lProcCount; ii++)
  {
    if (std::abs(lExpectedPerRank - (double) lNetTop->FluidSitesOnEachProcessor[ii]) > 1.0)
    {
      printf(
             "Domain decomposition problem: Expecting %f sites per processor, but rank %d had %d.\n",
             lExpectedPerRank, ii, lNetTop->FluidSitesOnEachProcessor[ii]);
    }
  }

  delete[] lName;
  delete lNetTop;
}

void Tests::MakeAGlobLatDat(hemelb::lb::GlobalLatticeData & oGlobLatDat)
{
  oGlobLatDat.SetBasicDetails(8, 8, 8, 4);

  for (unsigned int x = 0; x < oGlobLatDat.GetXBlockCount(); x++)
  {
    for (unsigned int y = 0; y < oGlobLatDat.GetYBlockCount(); y++)
    {
      for (unsigned int z = 0; z < oGlobLatDat.GetZBlockCount(); z++)
      {
        oGlobLatDat .Blocks[x * 64 + y * 8 + z].ProcessorRankForEachBlockSite
            = new int[oGlobLatDat.SitesPerBlockVolumeUnit];
        for (unsigned int lSite = 0; lSite < oGlobLatDat.SitesPerBlockVolumeUnit; lSite++)
        {
          oGlobLatDat .Blocks[x * 64 + y * 8 + z].ProcessorRankForEachBlockSite[lSite] = -1;
        }
      }
    }
  }
}
