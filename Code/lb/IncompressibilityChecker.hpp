// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_LB_INCOMPRESSIBILITYCHECKER_HPP
#define HEMELB_LB_INCOMPRESSIBILITYCHECKER_HPP

#include "lb/IncompressibilityChecker.h"

namespace hemelb
{
  namespace lb
  {

    template<class BroadcastPolicy>
    IncompressibilityChecker<BroadcastPolicy>::IncompressibilityChecker(const geometry::LatticeData * latticeData,
                                                                        net::Net* net,
                                                                        SimulationState* simState,
                                                                        lb::MacroscopicPropertyCache& propertyCache,
                                                                        reporting::Timers& timings,
                                                                        distribn_t maximumRelativeDensityDifferenceAllowed) :
        BroadcastPolicy(net, simState, SPREADFACTOR), mLatDat(latticeData), propertyCache(propertyCache), mSimState(simState), timings(timings), maximumRelativeDensityDifferenceAllowed(maximumRelativeDensityDifferenceAllowed), globalDensityTracker(NULL)
    {
      /*
       *  childrenDensitiesSerialised must be initialised to something sensible since ReceiveFromChildren won't
       *  fill it in completely unless the logarithm base SPREADFACTOR of the number of processes is an integer.
       */
      for (unsigned leaf_index = 0; leaf_index < SPREADFACTOR; leaf_index++)
      {
        unsigned offset = leaf_index * DensityTracker::DENSITY_TRACKER_SIZE;
        for (unsigned tracker_index = 0; tracker_index < DensityTracker::DENSITY_TRACKER_SIZE; tracker_index++)
        {
          switch (tracker_index)
          {
            case DensityTracker::MIN_DENSITY:
            case DensityTracker::MAX_DENSITY:
              childrenDensitiesSerialised[offset + tracker_index] = REFERENCE_DENSITY;
              break;
            case DensityTracker::MAX_VELOCITY_MAGNITUDE:
              childrenDensitiesSerialised[offset + tracker_index] = 0.0;
              break;
            default:
              // This should never trip. It only occurs when a new entry is added to the density tracker and
              // no suitable initialisation is provided.
              assert(false);
          }
        }
      }

    }

    template<class BroadcastPolicy>
    IncompressibilityChecker<BroadcastPolicy>::~IncompressibilityChecker()
    {
    }

    template<class BroadcastPolicy>
    distribn_t IncompressibilityChecker<BroadcastPolicy>::GetGlobalSmallestDensity() const
    {
      assert(AreDensitiesAvailable());
      return (*globalDensityTracker)[DensityTracker::MIN_DENSITY];
    }

    template<class BroadcastPolicy>
    distribn_t IncompressibilityChecker<BroadcastPolicy>::GetGlobalLargestDensity() const
    {
      assert(AreDensitiesAvailable());
      return (*globalDensityTracker)[DensityTracker::MAX_DENSITY];
    }

    template<class BroadcastPolicy>
    double IncompressibilityChecker<BroadcastPolicy>::GetMaxRelativeDensityDifference() const
    {
      distribn_t maxDensityDiff = GetGlobalLargestDensity() - GetGlobalSmallestDensity();
      assert(maxDensityDiff >= 0.0);
      return maxDensityDiff / REFERENCE_DENSITY;
    }

    template<class BroadcastPolicy>
    double IncompressibilityChecker<BroadcastPolicy>::GetMaxRelativeDensityDifferenceAllowed() const
    {
      return maximumRelativeDensityDifferenceAllowed;
    }

    template<class BroadcastPolicy>
    void IncompressibilityChecker<BroadcastPolicy>::PostReceiveFromChildren(unsigned long splayNumber)
    {
      timings[hemelb::reporting::Timers::monitoring].Start();

      for (int childIndex = 0; childIndex < (int) SPREADFACTOR; childIndex++)
      {
        DensityTracker childDensities(&childrenDensitiesSerialised[childIndex * DensityTracker::DENSITY_TRACKER_SIZE]);

        upwardsDensityTracker.UpdateDensityTracker(childDensities);
      }

      timings[hemelb::reporting::Timers::monitoring].Stop();
    }

    template<class BroadcastPolicy>
    void IncompressibilityChecker<BroadcastPolicy>::PostSendToParent(unsigned long splayNumber)
    {
      timings[hemelb::reporting::Timers::monitoring].Start();

      for (site_t i = 0; i < mLatDat->GetLocalFluidSiteCount(); i++)
      {
        upwardsDensityTracker.UpdateDensityTracker(propertyCache.densityCache.Get(i),
                                                   propertyCache.velocityCache.Get(i).GetMagnitude());
      }

      timings[hemelb::reporting::Timers::monitoring].Stop();
    }

    template<class BroadcastPolicy>
    void IncompressibilityChecker<BroadcastPolicy>::ProgressFromChildren(unsigned long splayNumber)
    {
      this->ReceiveFromChildren(childrenDensitiesSerialised, DensityTracker::DENSITY_TRACKER_SIZE);
    }

    template<class BroadcastPolicy>
    void IncompressibilityChecker<BroadcastPolicy>::ProgressFromParent(unsigned long splayNumber)
    {
      this->ReceiveFromParent(downwardsDensityTracker.GetDensitiesArray(), DensityTracker::DENSITY_TRACKER_SIZE);
    }

    template<class BroadcastPolicy>
    void IncompressibilityChecker<BroadcastPolicy>::ProgressToChildren(unsigned long splayNumber)
    {
      this->SendToChildren(downwardsDensityTracker.GetDensitiesArray(), DensityTracker::DENSITY_TRACKER_SIZE);
    }

    template<class BroadcastPolicy>
    void IncompressibilityChecker<BroadcastPolicy>::ProgressToParent(unsigned long splayNumber)
    {
      this->SendToParent(upwardsDensityTracker.GetDensitiesArray(), DensityTracker::DENSITY_TRACKER_SIZE);
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

    template<class BroadcastPolicy>
    bool IncompressibilityChecker<BroadcastPolicy>::AreDensitiesAvailable() const
    {
      return (globalDensityTracker != NULL);
    }

    template<class BroadcastPolicy>
    bool IncompressibilityChecker<BroadcastPolicy>::IsDensityDiffWithinRange() const
    {
      return (GetMaxRelativeDensityDifference() < maximumRelativeDensityDifferenceAllowed);
    }

    template<class BroadcastPolicy>
    void IncompressibilityChecker<BroadcastPolicy>::Report(ctemplate::TemplateDictionary& dictionary)
    {
      if (AreDensitiesAvailable() && !IsDensityDiffWithinRange())
      {
        ctemplate::TemplateDictionary *incomp = dictionary.AddSectionDictionary("DENSITIES");
        incomp->SetFormattedValue("ALLOWED", "%.1f%%", GetMaxRelativeDensityDifferenceAllowed() * 100);
        incomp->SetFormattedValue("ACTUAL", "%.1f%%", GetMaxRelativeDensityDifference() * 100);
      }
    }

    template<class BroadcastPolicy>
    double IncompressibilityChecker<BroadcastPolicy>::GetGlobalLargestVelocityMagnitude() const
    {
      assert(AreDensitiesAvailable());
      return (*globalDensityTracker)[DensityTracker::MAX_VELOCITY_MAGNITUDE];
    }
  }
}

#endif /* HEMELB_LB_INCOMPRESSIBILITYCHECKER_HPP */
