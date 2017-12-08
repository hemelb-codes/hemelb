
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_ADELB_HPP
#define HEMELB_LB_ADELB_HPP

#include "io/writers/xdr/XdrMemWriter.h"
#include "lb/ADELB.h"

namespace hemelb
{
  namespace lb
  {

    template<class LatticeType>
    hemelb::lb::LbmParameters* ADELBM<LatticeType>::GetLbmParams()
    {
      return &mParams;
    }

    template<class LatticeType>
    lb::MacroscopicPropertyCache& ADELBM<LatticeType>::GetPropertyCache()
    {
      return propertyCache;
    }

    template<class LatticeType>
    ADELBM<LatticeType>::ADELBM(configuration::SimConfig *iSimulationConfig,
                                net::Net* net,
                                geometry::LatticeData* latDat,
                                SimulationState* simState,
                                reporting::Timers &atimings,
                                geometry::neighbouring::NeighbouringDataManager *neighbouringDataManager,
                                lb::MacroscopicPropertyCache& coupledPropertyCache) :
      mSimConfig(iSimulationConfig), mNet(net), mLatDat(latDat), mState(simState), 
          mParams(iSimulationConfig->GetTimeStepLength(), iSimulationConfig->GetVoxelSize()), timings(atimings),
          propertyCache(*simState, *latDat), neighbouringDataManager(neighbouringDataManager), mCoupledPropertyCache(coupledPropertyCache)
    {
      ReadParameters();
    }

    template<class LatticeType>
    void ADELBM<LatticeType>::CalculateMouseFlowField(const ScreenDensity densityIn,
                                                      const LatticeDensity density_threshold_min,
                                                      const LatticeDensity density_threshold_minmax_inv,
                                                      PhysicalDensity &mouse_density)
    {
      LatticeDensity density = density_threshold_min + densityIn / density_threshold_minmax_inv;

      mouse_density = mUnits->ConvertDensityToPhysicalUnits(density);
    }

    template<class LatticeType>
    void ADELBM<LatticeType>::InitInitParamsSiteRanges(kernels::InitParams& initParams, unsigned& state)
    {
      initParams.siteRanges.resize(2);

      initParams.siteRanges[0].first = 0;
      initParams.siteRanges[1].first = mLatDat->GetMidDomainSiteCount();
      state = 0;
      initParams.siteRanges[0].second = initParams.siteRanges[0].first + mLatDat->GetMidDomainCollisionCount(state);
      initParams.siteRanges[1].second = initParams.siteRanges[1].first + mLatDat->GetDomainEdgeCollisionCount(state);

      initParams.siteCount = mLatDat->GetMidDomainCollisionCount(state) + mLatDat->GetDomainEdgeCollisionCount(state);
    }

    template<class LatticeType>
    void ADELBM<LatticeType>:: AdvanceInitParamsSiteRanges(kernels::InitParams& initParams, unsigned& state)
    {
      initParams.siteRanges[0].first += mLatDat->GetMidDomainCollisionCount(state);
      initParams.siteRanges[1].first += mLatDat->GetDomainEdgeCollisionCount(state);
      ++state;
      initParams.siteRanges[0].second = initParams.siteRanges[0].first + mLatDat->GetMidDomainCollisionCount(state);
      initParams.siteRanges[1].second = initParams.siteRanges[1].first + mLatDat->GetDomainEdgeCollisionCount(state);

      initParams.siteCount = mLatDat->GetMidDomainCollisionCount(state) + mLatDat->GetDomainEdgeCollisionCount(state);
    }

    template<class LatticeType>
    void ADELBM<LatticeType>::InitCollisions()
    {
      /**
       * Ensure the boundary objects have all info necessary.
       */
      PrepareBoundaryObjects();

      // TODO Note that the convergence checking is not yet implemented in the
      // new boundary condition hierarchy system.
      // It'd be nice to do this with something like
      // MidFluidCollision = new ConvergenceCheckingWrapper(new WhateverMidFluidCollision());

      kernels::InitParams initParams = kernels::InitParams();
      initParams.latDat = mLatDat;
      initParams.lbmParams = &mParams;
      initParams.neighbouringDataManager = neighbouringDataManager;

      unsigned collId;
      InitInitParamsSiteRanges(initParams, collId);
      mMidFluidCollision = new tMidFluidCollision(initParams);

      AdvanceInitParamsSiteRanges(initParams, collId);
      initParams.advectionDiffusionBoundaryObject = mStentValues;
      mWallCollision = new tWallCollision(initParams);

      AdvanceInitParamsSiteRanges(initParams, collId);
      initParams.boundaryObject = mInletValues;
      mInletCollision = new tInletCollision(initParams);

      AdvanceInitParamsSiteRanges(initParams, collId);
      initParams.boundaryObject = mOutletValues;
      mOutletCollision = new tOutletCollision(initParams);

      AdvanceInitParamsSiteRanges(initParams, collId);
      initParams.advectionDiffusionBoundaryObject = mStentValues;
      initParams.boundaryObject = mInletValues;
      mInletWallCollision = new tInletWallCollision(initParams);

      AdvanceInitParamsSiteRanges(initParams, collId);
      initParams.advectionDiffusionBoundaryObject = mStentValues;
      initParams.boundaryObject = mOutletValues;
      mOutletWallCollision = new tOutletWallCollision(initParams);
    }

    template<class LatticeType>
    void ADELBM<LatticeType>::Initialise(vis::Control* iControl,
                                         stents::BoundaryValues* iStentValues,
                                         iolets::BoundaryValues* iOutletValues,
                                         iolets::BoundaryValues* iInletValues,
                                         const util::UnitConverter* iUnits)
    {
      mStentValues = iStentValues;
      mOutletValues = iOutletValues;
      mInletValues = iInletValues;
      mUnits = iUnits;

      InitCollisions();

      SetInitialConditions();

      mVisControl = iControl;
    }

    template<class LatticeType>
    void ADELBM<LatticeType>::PrepareBoundaryObjects()
    {
      // First, iterate through all of the inlet and outlet objects, finding out the minimum density seen in the simulation.
      distribn_t minDensity = std::numeric_limits<distribn_t>::max();

      for (unsigned stent = 0; stent < mStentValues->GetLocalStentCount(); ++stent)
      {
        minDensity = std::min(minDensity, mStentValues->GetLocalStent(stent)->GetDensityMin());
      } 

      // Now go through them again, informing them of the minimum density.
      for (unsigned stent = 0; stent < mStentValues->GetLocalStentCount(); ++stent)
      {
        mStentValues->GetLocalStent(stent)->SetMinimumSimulationDensity(minDensity);
      }
    }

    template<class LatticeType>
    void ADELBM<LatticeType>::SetInitialConditions()
    {
      distribn_t density = mUnits->ConvertDensityToLatticeUnits(mSimConfig->GetInitialDensity());

      for (site_t i = 0; i < mLatDat->GetLocalFluidSiteCount(); i++)
      {
        distribn_t f_eq[LatticeType::NUMVECTORS];

        LatticeType::CalculateADEFeq(density, 0.0, 0.0, 0.0, f_eq);

        distribn_t* f_old_p = mLatDat->GetFOld(i * LatticeType::NUMVECTORS);
        distribn_t* f_new_p = mLatDat->GetFNew(i * LatticeType::NUMVECTORS);

        for (unsigned int l = 0; l < LatticeType::NUMVECTORS; l++)
        {
          f_new_p[l] = f_old_p[l] = f_eq[l];
        }
      }
    }

    template<class LatticeType>
    void ADELBM<LatticeType>::RequestComms()
    {
      timings[hemelb::reporting::Timers::lb].Start();

      // Delegate to the lattice data object to post the asynchronous sends and receives
      // (via the Net object).
      // NOTE that this doesn't actually *perform* the sends and receives, it asks the Net
      // to include them in the ISends and IRecvs that happen later.
      mLatDat->SendAndReceive(mNet);

      timings[hemelb::reporting::Timers::lb].Stop();
    }

    template<class LatticeType>
    void ADELBM<LatticeType>::PreSend()
    {
      timings[hemelb::reporting::Timers::lb].Start();
      timings[hemelb::reporting::Timers::lb_calc].Start();

      /**
       * In the PreSend phase, we do LB on all the sites that need to have results sent to
       * neighbouring ranks ('domainEdge' sites). In site id terms, this means we start at the
       * end of the sites whose neighbours all lie on this rank ('midDomain'), then progress
       * through the sites of each type in turn.
       */
      site_t offset = mLatDat->GetMidDomainSiteCount();

      StreamAndCollide(mMidFluidCollision, offset, mLatDat->GetDomainEdgeCollisionCount(0), mCoupledPropertyCache);
      offset += mLatDat->GetDomainEdgeCollisionCount(0);
      
      mStentValues->FinishReceive();
      StreamAndCollide(mWallCollision, offset, mLatDat->GetDomainEdgeCollisionCount(1), mCoupledPropertyCache);
      offset += mLatDat->GetDomainEdgeCollisionCount(1);

      mInletValues->FinishReceive();
      StreamAndCollide(mInletCollision, offset, mLatDat->GetDomainEdgeCollisionCount(2), mCoupledPropertyCache);
      offset += mLatDat->GetDomainEdgeCollisionCount(2);

      mOutletValues->FinishReceive();
      StreamAndCollide(mOutletCollision, offset, mLatDat->GetDomainEdgeCollisionCount(3), mCoupledPropertyCache);
      offset += mLatDat->GetDomainEdgeCollisionCount(3);

      StreamAndCollide(mInletWallCollision, offset, mLatDat->GetDomainEdgeCollisionCount(4), mCoupledPropertyCache);
      offset += mLatDat->GetDomainEdgeCollisionCount(4);

      StreamAndCollide(mOutletWallCollision, offset, mLatDat->GetDomainEdgeCollisionCount(5), mCoupledPropertyCache);

      timings[hemelb::reporting::Timers::lb_calc].Stop();
      timings[hemelb::reporting::Timers::lb].Stop();
    }

    template<class LatticeType>
    void ADELBM<LatticeType>::PreReceive()
    {
      timings[hemelb::reporting::Timers::lb].Start();
      timings[hemelb::reporting::Timers::lb_calc].Start();

      /**
       * In the PreReceive phase, we perform LB for all the sites whose neighbours lie on this
       * rank ('midDomain' rather than 'domainEdge' sites). Ideally this phase is the longest bit (maximising time for the asynchronous sends
       * and receives to complete).
       *
       * In site id terms, this means starting at the first site and progressing through the
       * midDomain sites, one type at a time.
       */
      site_t offset = 0;

      StreamAndCollide(mMidFluidCollision, offset, mLatDat->GetMidDomainCollisionCount(0), mCoupledPropertyCache);
      offset += mLatDat->GetMidDomainCollisionCount(0);

      StreamAndCollide(mWallCollision, offset, mLatDat->GetMidDomainCollisionCount(1), mCoupledPropertyCache);
      offset += mLatDat->GetMidDomainCollisionCount(1);

      StreamAndCollide(mInletCollision, offset, mLatDat->GetMidDomainCollisionCount(2), mCoupledPropertyCache);
      offset += mLatDat->GetMidDomainCollisionCount(2);

      StreamAndCollide(mOutletCollision, offset, mLatDat->GetMidDomainCollisionCount(3), mCoupledPropertyCache);
      offset += mLatDat->GetMidDomainCollisionCount(3);

      StreamAndCollide(mInletWallCollision, offset, mLatDat->GetMidDomainCollisionCount(4), mCoupledPropertyCache);
      offset += mLatDat->GetMidDomainCollisionCount(4);

      StreamAndCollide(mOutletWallCollision, offset, mLatDat->GetMidDomainCollisionCount(5), mCoupledPropertyCache);

      timings[hemelb::reporting::Timers::lb_calc].Stop();
      timings[hemelb::reporting::Timers::lb].Stop();
    }

    template<class LatticeType>
    void ADELBM<LatticeType>::PostReceive()
    {
      timings[hemelb::reporting::Timers::lb].Start();

      // Copy the distribution functions received from the neighbouring
      // processors into the destination buffer "f_new".
      // This is done here, after receiving the sent distributions from neighbours.
      mLatDat->CopyReceived();

      // Do any cleanup steps necessary on boundary nodes
      site_t offset = mLatDat->GetMidDomainSiteCount();

      timings[hemelb::reporting::Timers::lb_calc].Start();

      //TODO yup, this is horrible. If you read this, please improve the following code.
      PostStep(mMidFluidCollision, offset, mLatDat->GetDomainEdgeCollisionCount(0));
      offset += mLatDat->GetDomainEdgeCollisionCount(0);

      PostStep(mWallCollision, offset, mLatDat->GetDomainEdgeCollisionCount(1));
      offset += mLatDat->GetDomainEdgeCollisionCount(1);

      PostStep(mInletCollision, offset, mLatDat->GetDomainEdgeCollisionCount(2));
      offset += mLatDat->GetDomainEdgeCollisionCount(2);

      PostStep(mOutletCollision, offset, mLatDat->GetDomainEdgeCollisionCount(3));
      offset += mLatDat->GetDomainEdgeCollisionCount(3);

      PostStep(mInletWallCollision, offset, mLatDat->GetDomainEdgeCollisionCount(4));
      offset += mLatDat->GetDomainEdgeCollisionCount(4);

      PostStep(mOutletWallCollision, offset, mLatDat->GetDomainEdgeCollisionCount(5));

      offset = 0;

      PostStep(mMidFluidCollision, offset, mLatDat->GetMidDomainCollisionCount(0));
      offset += mLatDat->GetMidDomainCollisionCount(0);

      PostStep(mWallCollision, offset, mLatDat->GetMidDomainCollisionCount(1));
      offset += mLatDat->GetMidDomainCollisionCount(1);

      PostStep(mInletCollision, offset, mLatDat->GetMidDomainCollisionCount(2));
      offset += mLatDat->GetMidDomainCollisionCount(2);

      PostStep(mOutletCollision, offset, mLatDat->GetMidDomainCollisionCount(3));
      offset += mLatDat->GetMidDomainCollisionCount(3);

      PostStep(mInletWallCollision, offset, mLatDat->GetMidDomainCollisionCount(4));
      offset += mLatDat->GetMidDomainCollisionCount(4);

      PostStep(mOutletWallCollision, offset, mLatDat->GetMidDomainCollisionCount(5));

      timings[hemelb::reporting::Timers::lb_calc].Stop();
      timings[hemelb::reporting::Timers::lb].Stop();
    }

    template<class LatticeType>
    void ADELBM<LatticeType>::EndIteration()
    {
      timings[hemelb::reporting::Timers::lb].Start();
      timings[hemelb::reporting::Timers::lb_calc].Start();

      // Swap f_old and f_new ready for the next timestep.
      mLatDat->SwapOldAndNew();

      timings[hemelb::reporting::Timers::lb_calc].Stop();
      timings[hemelb::reporting::Timers::lb].Stop();
    }

    template<class LatticeType>
    ADELBM<LatticeType>::~ADELBM()
    {
      // Delete the collision and stream objects we've been using
      delete mMidFluidCollision;
      delete mWallCollision;
      delete mInletCollision;
      delete mOutletCollision;
      delete mInletWallCollision;
      delete mOutletWallCollision;
    }

    template<class LatticeType>
    void ADELBM<LatticeType>::ReadParameters()
    {
      std::vector<lb::stents::Stent*> stents = mSimConfig->GetStents();
      stentCount = stents.size();
    }

    template<class LatticeType>
    void ADELBM<LatticeType>::ReadVisParameters()
    {
      distribn_t density_min = std::numeric_limits<distribn_t>::max();
      distribn_t density_max = std::numeric_limits<distribn_t>::min();

      distribn_t velocity_max = mUnits->ConvertVelocityToLatticeUnits(mSimConfig->GetMaximumVelocity());

      for (int i = 0; i < StentCount(); i++)
      {
        density_min = util::NumericalFunctions::min(density_min, mStentValues->GetDensityMin(i));
        density_max = util::NumericalFunctions::max(density_max, mStentValues->GetDensityMax(i));
      }

      distribn_t lDensity_threshold_min = density_min;
      distribn_t lDensity_threshold_minmax_inv = 1.0F / (density_max - density_min);
      distribn_t lVelocity_threshold_max_inv = 1.0F / velocity_max;

      mVisControl->SetSomeParams(mSimConfig->GetVisualisationBrightness(),
                                 lDensity_threshold_min,
                                 lDensity_threshold_minmax_inv,
                                 lVelocity_threshold_max_inv,
                                 1.0);
    }
  }
}

#endif /* HEMELB_LB_ADELB_HPP */