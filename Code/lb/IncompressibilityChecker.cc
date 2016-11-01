// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "lb/IncompressibilityChecker.h"

namespace hemelb
{
  namespace comm
  {
    template<>
    MPI_Datatype MpiDataTypeTraits<lb::DensityEtc>::RegisterMpiDataType()
    {
      MPI_Datatype dType = lb::IncompressibilityChecker::GetDensityMpiDatatype();
      HEMELB_MPI_CALL(MPI_Type_commit, (&dType));
      return dType;
    }
  }

  namespace lb
  {
    /**
     * Create the MPI Datatype for sending DensityEtc structs.
     */
    MPI_Datatype IncompressibilityChecker::GetDensityMpiDatatype()
    {
      HEMELB_MPI_TYPE_BEGIN(dType, ::hemelb::lb::DensityEtc, 3);
      HEMELB_MPI_TYPE_ADD_MEMBER(min);
      HEMELB_MPI_TYPE_ADD_MEMBER(max);
      HEMELB_MPI_TYPE_ADD_MEMBER(maxVel);
      HEMELB_MPI_TYPE_END(dType, ::hemelb::lb::DensityEtc);
      return dType;
    }

    void IncompressibilityChecker::UpdateDensityEtc(const DensityEtc& in, DensityEtc& inout)
    {
      inout.min = std::min(in.min, inout.min);
      inout.max = std::max(in.max, inout.max);
      inout.maxVel = std::max(in.maxVel, inout.maxVel);
    }

    void IncompressibilityChecker::MpiOpUpdateFunc(void* invec, void* inoutvec, int *len,
                                                   MPI_Datatype *datatype)
    {
      const DensityEtc* iv = static_cast<const DensityEtc*>(invec);
      DensityEtc* iov = static_cast<DensityEtc*>(inoutvec);
      for (int i = 0; i < *len; ++i)
      {
        UpdateDensityEtc(iv[i], iov[i]);
      }
    }

    IncompressibilityChecker::IncompressibilityChecker(const geometry::LatticeData * latticeData,
						       comm::Communicator::ConstPtr comms,
						       SimulationState* simState,
						       lb::MacroscopicPropertyCache& propertyCache,
						       reporting::Timers& timings,
						       distribn_t maximumRelativeDensityDifferenceAllowed) :
        timestep::CollectiveActor(comms, timings[reporting::Timers::mpiWait]),
            mLatDat(latticeData), propertyCache(propertyCache), mSimState(simState),
            maximumRelativeDensityDifferenceAllowed(maximumRelativeDensityDifferenceAllowed),
            workTimer(timings[reporting::Timers::monitoring])
    {
      HEMELB_MPI_CALL(MPI_Op_create, (&IncompressibilityChecker::MpiOpUpdateFunc, 1, &reduction));
      localDensity.min = std::numeric_limits<double>::max();
      localDensity.max = -std::numeric_limits<double>::max();
      localDensity.maxVel = 0;
      globalDensity = localDensity;
    }

    IncompressibilityChecker::~IncompressibilityChecker()
    {
      HEMELB_MPI_CALL(MPI_Op_free, (&reduction));
    }

    distribn_t IncompressibilityChecker::GetGlobalSmallestDensity() const
    {
      return globalDensity.min;
    }

    distribn_t IncompressibilityChecker::GetGlobalLargestDensity() const
    {
      return globalDensity.max;
    }

    double IncompressibilityChecker::GetMaxRelativeDensityDifference() const
    {
      distribn_t maxDensityDiff = GetGlobalLargestDensity() - GetGlobalSmallestDensity();
      if (maxDensityDiff < 0.0)
        return std::numeric_limits<double>::quiet_NaN();

      return maxDensityDiff / REFERENCE_DENSITY;
    }

    double IncompressibilityChecker::GetMaxRelativeDensityDifferenceAllowed() const
    {
      return maximumRelativeDensityDifferenceAllowed;
    }

    void IncompressibilityChecker::PreSend()
    {
      for (site_t i = 0; i < mLatDat->GetLocalFluidSiteCount(); i++)
      {
        DensityEtc etc;
        etc.min = propertyCache.densityCache.Get(i);
        etc.max = propertyCache.densityCache.Get(i);
        etc.maxVel = propertyCache.velocityCache.Get(i).GetMagnitude();
        UpdateDensityEtc(etc, localDensity);
      }

    }
    /**
     * Initiate the collective.
     */
    void IncompressibilityChecker::Send(void)
    {
      // Begin collective.
      collectiveReq = collectiveComm->Iallreduce(localDensity, reduction, globalDensity);
    }

    bool IncompressibilityChecker::IsDensityDiffWithinRange() const
    {
      return (GetMaxRelativeDensityDifference() < maximumRelativeDensityDifferenceAllowed);
    }

    void IncompressibilityChecker::Report(ctemplate::TemplateDictionary& dictionary)
    {
      if (!IsDensityDiffWithinRange())
      {
        ctemplate::TemplateDictionary *incomp = dictionary.AddSectionDictionary("DENSITIES");
        incomp->SetFormattedValue("ALLOWED",
                                  "%.1f%%",
                                  GetMaxRelativeDensityDifferenceAllowed() * 100);
        incomp->SetFormattedValue("ACTUAL", "%.1f%%", GetMaxRelativeDensityDifference() * 100);
      }
    }

    double IncompressibilityChecker::GetGlobalLargestVelocityMagnitude() const
    {
      return globalDensity.maxVel;
    }

  }
}
