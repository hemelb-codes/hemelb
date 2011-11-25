#include "lb/IncompressibilityChecker.h"
#include "net/BroadcastMock.h"

namespace hemelb
{
  namespace lb
  {
    template<class BroadcastPolicy>
    IncompressibilityChecker<BroadcastPolicy>::DensityTracker::DensityTracker() :
      allocatedHere(true)
    {
      densitiesArray = new distribn_t[DENSITY_TRACKER_SIZE];

      /// TODO: do we have a policy on floating point constants?
      densitiesArray[MIN_DENSITY] = DBL_MAX;
      densitiesArray[MAX_DENSITY] = -DBL_MAX;
    }

    template<class BroadcastPolicy>
    IncompressibilityChecker<BroadcastPolicy>::DensityTracker::DensityTracker(distribn_t* const densityValues) :
      densitiesArray(densityValues), allocatedHere(false)
    {
    }

    template<class BroadcastPolicy>
    IncompressibilityChecker<BroadcastPolicy>::DensityTracker::~DensityTracker()
    {
      if (allocatedHere)
      {
        delete[] densitiesArray;
      }
    }

    template<class BroadcastPolicy>
    void IncompressibilityChecker<BroadcastPolicy>::DensityTracker::operator=(const DensityTracker& newValues)
    {
      for (unsigned trackerEntry = 0; trackerEntry < DENSITY_TRACKER_SIZE; trackerEntry++)
      {
        densitiesArray[trackerEntry] = newValues.GetDensitiesArray()[trackerEntry];
      }
    }

    template<class BroadcastPolicy>
    distribn_t& IncompressibilityChecker<BroadcastPolicy>::DensityTracker::operator[](DensityTrackerIndices densityIndex) const
    {
      return densitiesArray[densityIndex];
    }

    template<class BroadcastPolicy>
    distribn_t* IncompressibilityChecker<BroadcastPolicy>::DensityTracker::GetDensitiesArray() const
    {
      return densitiesArray;
    }

    template<class BroadcastPolicy>
    void IncompressibilityChecker<BroadcastPolicy>::DensityTracker::UpdateDensityTracker(const DensityTracker& newValues)
    {
      if (newValues[MIN_DENSITY] < densitiesArray[MIN_DENSITY])
      {
        densitiesArray[MIN_DENSITY] = newValues[MIN_DENSITY];
      }
      if (newValues[MAX_DENSITY] > densitiesArray[MAX_DENSITY])
      {
        densitiesArray[MAX_DENSITY] = newValues[MAX_DENSITY];
      }
    }

    template<class BroadcastPolicy>
    void IncompressibilityChecker<BroadcastPolicy>::DensityTracker::UpdateDensityTracker(distribn_t newValue)
    {
      if (newValue < densitiesArray[MIN_DENSITY])
      {
        densitiesArray[MIN_DENSITY] = newValue;
      }
      if (newValue > densitiesArray[MAX_DENSITY])
      {
        densitiesArray[MAX_DENSITY] = newValue;
      }
    }

    template<class BroadcastPolicy>
    IncompressibilityChecker<BroadcastPolicy>::IncompressibilityChecker(const geometry::LatticeData * latticeData,
                                                                        net::Net* net,
                                                                        SimulationState* simState) :
      BroadcastPolicy(net, simState, SPREADFACTOR), mLatDat(latticeData), mSimState(simState),
          globalDensityTracker(NULL)
    {
    }

    template<class BroadcastPolicy>
    IncompressibilityChecker<BroadcastPolicy>::~IncompressibilityChecker()
    {
    }

    template<class BroadcastPolicy>
    distribn_t IncompressibilityChecker<BroadcastPolicy>::GetGlobalSmallestDensity()
    {
      assert(globalDensityTracker != NULL);
      return (*globalDensityTracker)[DensityTracker::MIN_DENSITY];
    }

    template<class BroadcastPolicy>
    distribn_t IncompressibilityChecker<BroadcastPolicy>::GetGlobalLargestDensity()
    {
      assert(globalDensityTracker != NULL);
      return (*globalDensityTracker)[DensityTracker::MAX_DENSITY];
    }

    template<class BroadcastPolicy>
    distribn_t IncompressibilityChecker<BroadcastPolicy>::GetMaxDensityDifference()
    {
      distribn_t maxDiff = GetGlobalLargestDensity() - GetGlobalSmallestDensity();
      assert(maxDiff >= 0.0);
      return maxDiff;
    }

    template<class BroadcastPolicy>
    void IncompressibilityChecker<BroadcastPolicy>::PostReceiveFromChildren(unsigned long splayNumber)
    {
      for (int childIndex = 0; childIndex < (int) SPREADFACTOR; childIndex++)
      {
        DensityTracker childDensities(&childrenDensityTracker[childIndex
            * DensityTracker::DENSITY_TRACKER_SIZE]);

        upwardsDensityTracker.UpdateDensityTracker(childDensities);
      }
    }

    template<class BroadcastPolicy>
    void IncompressibilityChecker<BroadcastPolicy>::ProgressFromChildren(unsigned long splayNumber)
    {
      this->ReceiveFromChildren(&childrenDensityTracker[0], DensityTracker::DENSITY_TRACKER_SIZE);
    }

    template<class BroadcastPolicy>
    void IncompressibilityChecker<BroadcastPolicy>::ProgressFromParent(unsigned long splayNumber)
    {
      this->ReceiveFromParent(downwardsDensityTracker.GetDensitiesArray(),
                              DensityTracker::DENSITY_TRACKER_SIZE);
    }

    template<class BroadcastPolicy>
    void IncompressibilityChecker<BroadcastPolicy>::ProgressToChildren(unsigned long splayNumber)
    {
      this->SendToChildren(downwardsDensityTracker.GetDensitiesArray(),
                           DensityTracker::DENSITY_TRACKER_SIZE);
    }

    template<class BroadcastPolicy>
    void IncompressibilityChecker<BroadcastPolicy>::ProgressToParent(unsigned long splayNumber)
    {
      distribn_t localDensity;
      for (site_t i = 0; i < mLatDat->GetLocalFluidSiteCount(); i++)
      {
        /// TODO Refactor into a method in the lattice class that computes *just* the density
        localDensity = 0.0;
        for (unsigned int l = 0; l < D3Q15::NUMVECTORS; l++)
        {
          localDensity += *mLatDat->GetFNew(i * D3Q15::NUMVECTORS + l);
        }
        // End of refactor

        upwardsDensityTracker.UpdateDensityTracker(localDensity);
      }

      this->SendToParent(upwardsDensityTracker.GetDensitiesArray(),
                         DensityTracker::DENSITY_TRACKER_SIZE);
    }

    template<class BroadcastPolicy>
    void IncompressibilityChecker<BroadcastPolicy>::TopNodeAction()
    {
      downwardsDensityTracker = upwardsDensityTracker;
    }

    template<class BroadcastPolicy>
    void IncompressibilityChecker<BroadcastPolicy>::Effect()
    {
      globalDensityTracker = &downwardsDensityTracker;
    }

    // Explicit instantiation
    template class IncompressibilityChecker<net::PhasedBroadcastRegular<> > ;
    template class IncompressibilityChecker<net::BroadcastMock> ;

  }
}
