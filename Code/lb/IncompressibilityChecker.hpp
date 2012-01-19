#ifndef HEMELB_LB_INCOMPRESSIBILITYCHECKER_HPP
#define HEMELB_LB_INCOMPRESSIBILITYCHECKER_HPP

#include "lb/IncompressibilityChecker.h"

namespace hemelb
{
  namespace lb
  {
    template<class BroadcastPolicy, class Lattice>
    IncompressibilityChecker<BroadcastPolicy, Lattice>::DensityTracker::DensityTracker() :
        allocatedHere(true)
    {
      densitiesArray = new distribn_t[DENSITY_TRACKER_SIZE];

      //! @todo #23 do we have a policy on floating point constants?
      densitiesArray[MIN_DENSITY] = DBL_MAX;
      densitiesArray[MAX_DENSITY] = -DBL_MAX;
    }

    template<class BroadcastPolicy, class Lattice>
    IncompressibilityChecker<BroadcastPolicy, Lattice>::DensityTracker::DensityTracker(distribn_t* const densityValues) :
        densitiesArray(densityValues), allocatedHere(false)
    {
    }

    template<class BroadcastPolicy, class Lattice>
    IncompressibilityChecker<BroadcastPolicy, Lattice>::DensityTracker::~DensityTracker()
    {
      if (allocatedHere)
      {
        delete[] densitiesArray;
      }
    }

    template<class BroadcastPolicy, class Lattice>
    void IncompressibilityChecker<BroadcastPolicy, Lattice>::DensityTracker::operator=(const DensityTracker& newValues)
    {
      for (unsigned trackerEntry = 0; trackerEntry < DENSITY_TRACKER_SIZE; trackerEntry++)
      {
        densitiesArray[trackerEntry] = newValues.GetDensitiesArray()[trackerEntry];
      }
    }

    template<class BroadcastPolicy, class Lattice>
    distribn_t& IncompressibilityChecker<BroadcastPolicy, Lattice>::DensityTracker::operator[](DensityTrackerIndices densityIndex) const
    {
      return densitiesArray[densityIndex];
    }

    template<class BroadcastPolicy, class Lattice>
    distribn_t* IncompressibilityChecker<BroadcastPolicy, Lattice>::DensityTracker::GetDensitiesArray() const
    {
      return densitiesArray;
    }

    template<class BroadcastPolicy, class Lattice>
    void IncompressibilityChecker<BroadcastPolicy, Lattice>::DensityTracker::UpdateDensityTracker(const DensityTracker& newValues)
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

    template<class BroadcastPolicy, class Lattice>
    void IncompressibilityChecker<BroadcastPolicy, Lattice>::DensityTracker::UpdateDensityTracker(distribn_t newValue)
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

    template<class BroadcastPolicy, class Lattice>
    IncompressibilityChecker<BroadcastPolicy, Lattice>::IncompressibilityChecker(const geometry::LatticeData * latticeData,
                                                                                 net::Net* net,
                                                                                 SimulationState* simState,
                                                                                 reporting::Timers& timings,
                                                                                 distribn_t maximumRelativeDensityDifferenceAllowed) :
        BroadcastPolicy(net, simState, SPREADFACTOR), mLatDat(latticeData), mSimState(simState), timings(timings), maximumRelativeDensityDifferenceAllowed(maximumRelativeDensityDifferenceAllowed), globalDensityTracker(NULL)
    {
      /*
       *  childrenDensitiesSerialised must be initialised to something sensible since ReceiveFromChildren won't
       *  fill it in completely unless the logarithm base SPREADFACTOR of the number of processes is an integer.
       */
      for (unsigned int index = 0; index < SPREADFACTOR * DensityTracker::DENSITY_TRACKER_SIZE; index++)
      {
        childrenDensitiesSerialised[index] = REFERENCE_DENSITY;
      }

    }

    template<class BroadcastPolicy, class Lattice>
    IncompressibilityChecker<BroadcastPolicy, Lattice>::~IncompressibilityChecker()
    {
    }

    template<class BroadcastPolicy, class Lattice>
    distribn_t IncompressibilityChecker<BroadcastPolicy, Lattice>::GetGlobalSmallestDensity() const
    {
      assert(AreDensitiesAvailable());
      return (*globalDensityTracker)[DensityTracker::MIN_DENSITY];
    }

    template<class BroadcastPolicy, class Lattice>
    distribn_t IncompressibilityChecker<BroadcastPolicy, Lattice>::GetGlobalLargestDensity() const
    {
      assert(AreDensitiesAvailable());
      return (*globalDensityTracker)[DensityTracker::MAX_DENSITY];
    }

    template<class BroadcastPolicy, class Lattice>
    double IncompressibilityChecker<BroadcastPolicy, Lattice>::GetMaxRelativeDensityDifference() const
    {
      distribn_t maxDensityDiff = GetGlobalLargestDensity() - GetGlobalSmallestDensity();
      assert(maxDensityDiff >= 0.0);
      return maxDensityDiff / REFERENCE_DENSITY;
    }

    template<class BroadcastPolicy, class Lattice>
    double IncompressibilityChecker<BroadcastPolicy, Lattice>::GetMaxRelativeDensityDifferenceAllowed() const
    {
      return maximumRelativeDensityDifferenceAllowed;
    }

    template<class BroadcastPolicy, class Lattice>
    void IncompressibilityChecker<BroadcastPolicy, Lattice>::PostReceiveFromChildren(unsigned long splayNumber)
    {
      timings[hemelb::reporting::Timers::monitoring].Start();

      for (int childIndex = 0; childIndex < (int) SPREADFACTOR; childIndex++)
      {
        DensityTracker childDensities(&childrenDensitiesSerialised[childIndex * DensityTracker::DENSITY_TRACKER_SIZE]);

        upwardsDensityTracker.UpdateDensityTracker(childDensities);
      }

      timings[hemelb::reporting::Timers::monitoring].Stop();
    }

    template<class BroadcastPolicy, class Lattice>
    void IncompressibilityChecker<BroadcastPolicy, Lattice>::ProgressFromChildren(unsigned long splayNumber)
    {
      this->ReceiveFromChildren(childrenDensitiesSerialised, DensityTracker::DENSITY_TRACKER_SIZE);
    }

    template<class BroadcastPolicy, class Lattice>
    void IncompressibilityChecker<BroadcastPolicy, Lattice>::ProgressFromParent(unsigned long splayNumber)
    {
      this->ReceiveFromParent(downwardsDensityTracker.GetDensitiesArray(), DensityTracker::DENSITY_TRACKER_SIZE);
    }

    template<class BroadcastPolicy, class Lattice>
    void IncompressibilityChecker<BroadcastPolicy, Lattice>::ProgressToChildren(unsigned long splayNumber)
    {
      this->SendToChildren(downwardsDensityTracker.GetDensitiesArray(), DensityTracker::DENSITY_TRACKER_SIZE);
    }

    template<class BroadcastPolicy, class Lattice>
    void IncompressibilityChecker<BroadcastPolicy, Lattice>::ProgressToParent(unsigned long splayNumber)
    {
      timings[hemelb::reporting::Timers::monitoring].Start();

      distribn_t localDensity;
      for (site_t i = 0; i < mLatDat->GetLocalFluidSiteCount(); i++)
      {
        //! @todo #23 Refactor into a method in the lattice class that computes *just* the density
        localDensity = 0.0;
        for (unsigned int l = 0; l < Lattice::NUMVECTORS; l++)
        {
          localDensity += *mLatDat->GetFNew(i * Lattice::NUMVECTORS + l);
        }
        // End of refactor

        upwardsDensityTracker.UpdateDensityTracker(localDensity);
      }
      timings[hemelb::reporting::Timers::monitoring].Stop();

      this->SendToParent(upwardsDensityTracker.GetDensitiesArray(), DensityTracker::DENSITY_TRACKER_SIZE);
    }

    template<class BroadcastPolicy, class Lattice>
    void IncompressibilityChecker<BroadcastPolicy, Lattice>::TopNodeAction()
    {
      downwardsDensityTracker = upwardsDensityTracker;
    }

    template<class BroadcastPolicy, class Lattice>
    void IncompressibilityChecker<BroadcastPolicy, Lattice>::Effect()
    {
      globalDensityTracker = &downwardsDensityTracker;
    }

    template<class BroadcastPolicy, class Lattice>
    bool IncompressibilityChecker<BroadcastPolicy, Lattice>::AreDensitiesAvailable() const
    {
      return (globalDensityTracker != NULL);
    }

    template<class BroadcastPolicy, class Lattice>
    bool IncompressibilityChecker<BroadcastPolicy, Lattice>::IsDensityDiffWithinRange() const
    {
      return (GetMaxRelativeDensityDifference() < maximumRelativeDensityDifferenceAllowed);
    }

    template<class BroadcastPolicy, class Lattice>
    void IncompressibilityChecker<BroadcastPolicy, Lattice>::Report(ctemplate::TemplateDictionary& dictionary)
    {
      if (AreDensitiesAvailable() && !IsDensityDiffWithinRange())
      {
        ctemplate::TemplateDictionary *incomp = dictionary.AddSectionDictionary("DENSITIES");
        incomp->SetFormattedValue("ALLOWED", "%.1f%%", GetMaxRelativeDensityDifferenceAllowed() * 100);
        incomp->SetFormattedValue("ACTUAL", "%.1f%%", GetMaxRelativeDensityDifference() * 100);
      }
    }
  }
}

#endif /* HEMELB_LB_INCOMPRESSIBILITYCHECKER_HPP */
