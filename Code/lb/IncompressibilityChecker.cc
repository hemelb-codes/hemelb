// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#include "lb/lattices/D3Q15.h"
#include "lb/IncompressibilityChecker.hpp"

namespace hemelb
{
  namespace lb
  {
    DensityTracker::DensityTracker() :
        allocatedHere(true)
    {
      densitiesArray = new distribn_t[DENSITY_TRACKER_SIZE];

      //! @todo #23 do we have a policy on floating point constants?
      densitiesArray[MIN_DENSITY] = DBL_MAX;
      densitiesArray[MAX_DENSITY] = -DBL_MAX;
      densitiesArray[MAX_VELOCITY_MAGNITUDE] = 0.0;
    }

    DensityTracker::DensityTracker(distribn_t* const densityValues) :
        densitiesArray(densityValues), allocatedHere(false)
    {
    }

    DensityTracker::~DensityTracker()
    {
      if (allocatedHere)
      {
        delete[] densitiesArray;
      }
    }

    void DensityTracker::operator=(const DensityTracker& newValues)
    {
      for (unsigned trackerEntry = 0; trackerEntry < DENSITY_TRACKER_SIZE; trackerEntry++)
      {
        densitiesArray[trackerEntry] = newValues.GetDensitiesArray()[trackerEntry];
      }
    }

    distribn_t& DensityTracker::operator[](DensityTrackerIndices densityIndex) const
    {
      return densitiesArray[densityIndex];
    }

    distribn_t* DensityTracker::GetDensitiesArray() const
    {
      return densitiesArray;
    }

    void DensityTracker::UpdateDensityTracker(const DensityTracker& newValues)
    {
      if (newValues[MIN_DENSITY] < densitiesArray[MIN_DENSITY])
      {
        densitiesArray[MIN_DENSITY] = newValues[MIN_DENSITY];
      }
      if (newValues[MAX_DENSITY] > densitiesArray[MAX_DENSITY])
      {
        densitiesArray[MAX_DENSITY] = newValues[MAX_DENSITY];
      }
      if (newValues[MAX_VELOCITY_MAGNITUDE] > densitiesArray[MAX_VELOCITY_MAGNITUDE])
      {
        densitiesArray[MAX_VELOCITY_MAGNITUDE] = newValues[MAX_VELOCITY_MAGNITUDE];
      }
    }

    void DensityTracker::UpdateDensityTracker(distribn_t newDensity,
                                              distribn_t newVelocityMagnitude)
    {
      if (newDensity < densitiesArray[MIN_DENSITY])
      {
        densitiesArray[MIN_DENSITY] = newDensity;
      }
      if (newDensity > densitiesArray[MAX_DENSITY])
      {
        densitiesArray[MAX_DENSITY] = newDensity;
      }
      if (newVelocityMagnitude > densitiesArray[MAX_VELOCITY_MAGNITUDE])
      {
        densitiesArray[MAX_VELOCITY_MAGNITUDE] = newVelocityMagnitude;
      }
    }

    // Explicit instantiation
    template class IncompressibilityChecker<net::PhasedBroadcastRegular<> > ;
  }
}
