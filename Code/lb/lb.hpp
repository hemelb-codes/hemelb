// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_LB_HPP
#define HEMELB_LB_LB_HPP

#include "lb/lb.h"
#include "lb/InitialCondition.h"
#include "lb/InitialCondition.hpp"

namespace hemelb
{
  namespace lb
  {

    template<class TRAITS>
    hemelb::lb::LbmParameters* LBM<TRAITS>::GetLbmParams()
    {
      return &mParams;
    }

    template<class TRAITS>
    lb::MacroscopicPropertyCache& LBM<TRAITS>::GetPropertyCache()
    {
      return propertyCache;
    }

    template<class TRAITS>
    LBM<TRAITS>::LBM(LbmParameters iParams, net::Net* net,
                     geometry::FieldData* latDat, SimulationState* simState,
                     reporting::Timers &atimings,
                     geometry::neighbouring::NeighbouringDataManager *neighbouringDataManager) :
        mNet(net), mLatDat(latDat), mState(simState),
            mParams(std::move(iParams)),
            timings(atimings), propertyCache(*simState, latDat->GetDomain()),
            neighbouringDataManager(neighbouringDataManager)
    {
    }

    template<class TRAITS>
    void LBM<TRAITS>::InitInitParamsSiteRanges(kernels::InitParams& initParams, unsigned& state)
    {
      auto& dom = mLatDat->GetDomain();
      initParams.siteRanges.resize(2);

      initParams.siteRanges[0].first = 0;
      initParams.siteRanges[1].first = dom.GetMidDomainSiteCount();
      state = 0;
      initParams.siteRanges[0].second = initParams.siteRanges[0].first
          + dom.GetMidDomainCollisionCount(state);
      initParams.siteRanges[1].second = initParams.siteRanges[1].first
          + dom.GetDomainEdgeCollisionCount(state);

      initParams.siteCount = dom.GetMidDomainCollisionCount(state)
          + dom.GetDomainEdgeCollisionCount(state);
    }

    template<class TRAITS>
    void LBM<TRAITS>::AdvanceInitParamsSiteRanges(kernels::InitParams& initParams, unsigned& state)
    {
      auto& dom = mLatDat->GetDomain();
      initParams.siteRanges[0].first += dom.GetMidDomainCollisionCount(state);
      initParams.siteRanges[1].first += dom.GetDomainEdgeCollisionCount(state);
      ++state;
      initParams.siteRanges[0].second = initParams.siteRanges[0].first
          + dom.GetMidDomainCollisionCount(state);
      initParams.siteRanges[1].second = initParams.siteRanges[1].first
          + dom.GetDomainEdgeCollisionCount(state);

      initParams.siteCount = dom.GetMidDomainCollisionCount(state)
          + dom.GetDomainEdgeCollisionCount(state);
    }

    template<class TRAITS>
    void LBM<TRAITS>::InitCollisions()
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
      initParams.latDat = &mLatDat->GetDomain();
      initParams.lbmParams = &mParams;
      initParams.neighbouringDataManager = neighbouringDataManager;

      unsigned collId;
      InitInitParamsSiteRanges(initParams, collId);
      mMidFluidCollision = new tMidFluidCollision(initParams);

      AdvanceInitParamsSiteRanges(initParams, collId);
      mWallCollision = new tWallCollision(initParams);

      AdvanceInitParamsSiteRanges(initParams, collId);
      initParams.boundaryObject = mInletValues;
      mInletCollision = new tInletCollision(initParams);

      AdvanceInitParamsSiteRanges(initParams, collId);
      initParams.boundaryObject = mOutletValues;
      mOutletCollision = new tOutletCollision(initParams);

      AdvanceInitParamsSiteRanges(initParams, collId);
      initParams.boundaryObject = mInletValues;
      mInletWallCollision = new tInletWallCollision(initParams);

      AdvanceInitParamsSiteRanges(initParams, collId);
      initParams.boundaryObject = mOutletValues;
      mOutletWallCollision = new tOutletWallCollision(initParams);
    }

    template<class TRAITS>
    void LBM<TRAITS>::Initialise(iolets::BoundaryValues* iInletValues,
                                 iolets::BoundaryValues* iOutletValues)
    {
      mInletValues = iInletValues;
      mOutletValues = iOutletValues;

      InitCollisions();
    }

    template<class TRAITS>
    void LBM<TRAITS>::PrepareBoundaryObjects()
    {
      // First, iterate through all of the inlet and outlet objects, finding out the minimum density seen in the simulation.
      distribn_t minDensity = std::numeric_limits<distribn_t>::max();

      for (unsigned inlet = 0; inlet < mInletValues->GetLocalIoletCount(); ++inlet)
      {
        minDensity = std::min(minDensity, mInletValues->GetLocalIolet(inlet)->GetDensityMin());
      }

      for (unsigned outlet = 0; outlet < mOutletValues->GetLocalIoletCount(); ++outlet)
      {
        minDensity = std::min(minDensity, mOutletValues->GetLocalIolet(outlet)->GetDensityMin());
      }

      // Now go through them again, informing them of the minimum density.
      for (unsigned inlet = 0; inlet < mInletValues->GetLocalIoletCount(); ++inlet)
      {
        mInletValues->GetLocalIolet(inlet)->SetMinimumSimulationDensity(minDensity);
      }

      for (unsigned outlet = 0; outlet < mOutletValues->GetLocalIoletCount(); ++outlet)
      {
        mOutletValues->GetLocalIolet(outlet)->SetMinimumSimulationDensity(minDensity);
      }
    }

    template<class TRAITS>
    void LBM<TRAITS>::SetInitialConditions(lb::InitialCondition const& icond, const net::IOCommunicator& ioComms)
    {
      icond.SetFs<LatticeType>(mLatDat, ioComms);
      icond.SetTime(mState);
    }

    template<class TRAITS>
    void LBM<TRAITS>::RequestComms()
    {
      timings[hemelb::reporting::Timers::lb].Start();

      // Delegate to the lattice data object to post the asynchronous sends and receives
      // (via the Net object).
      // NOTE that this doesn't actually *perform* the sends and receives, it asks the Net
      // to include them in the ISends and IRecvs that happen later.
      mLatDat->SendAndReceive(mNet);

      timings[hemelb::reporting::Timers::lb].Stop();
    }

    template<class TRAITS>
    void LBM<TRAITS>::PreSend()
    {
      timings[hemelb::reporting::Timers::lb].Start();
      timings[hemelb::reporting::Timers::lb_calc].Start();

      /**
       * In the PreSend phase, we do LB on all the sites that need to have results sent to
       * neighbouring ranks ('domainEdge' sites). In site id terms, this means we start at the
       * end of the sites whose neighbours all lie on this rank ('midDomain'), then progress
       * through the sites of each type in turn.
       */
      auto& dom = mLatDat->GetDomain();
      site_t offset = dom.GetMidDomainSiteCount();

      log::Logger::Log<log::Debug, log::OnePerCore>("LBM - PreSend - StreamAndCollide");
      StreamAndCollide(mMidFluidCollision, offset, dom.GetDomainEdgeCollisionCount(0));
      offset += dom.GetDomainEdgeCollisionCount(0);

      StreamAndCollide(mWallCollision, offset, dom.GetDomainEdgeCollisionCount(1));
      offset += dom.GetDomainEdgeCollisionCount(1);

      mInletValues->FinishReceive();
      StreamAndCollide(mInletCollision, offset, dom.GetDomainEdgeCollisionCount(2));
      offset += dom.GetDomainEdgeCollisionCount(2);

      mOutletValues->FinishReceive();
      StreamAndCollide(mOutletCollision, offset, dom.GetDomainEdgeCollisionCount(3));
      offset += dom.GetDomainEdgeCollisionCount(3);

      StreamAndCollide(mInletWallCollision, offset, dom.GetDomainEdgeCollisionCount(4));
      offset += dom.GetDomainEdgeCollisionCount(4);

      StreamAndCollide(mOutletWallCollision, offset, dom.GetDomainEdgeCollisionCount(5));

      timings[hemelb::reporting::Timers::lb_calc].Stop();
      timings[hemelb::reporting::Timers::lb].Stop();
    }

    template<class TRAITS>
    void LBM<TRAITS>::PreReceive()
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
      auto& dom = mLatDat->GetDomain();
      site_t offset = 0;

      log::Logger::Log<log::Debug, log::OnePerCore>("LBM - PreReceive - StreamAndCollide");
      StreamAndCollide(mMidFluidCollision, offset, dom.GetMidDomainCollisionCount(0));
      offset += dom.GetMidDomainCollisionCount(0);

      StreamAndCollide(mWallCollision, offset, dom.GetMidDomainCollisionCount(1));
      offset += dom.GetMidDomainCollisionCount(1);

      StreamAndCollide(mInletCollision, offset, dom.GetMidDomainCollisionCount(2));
      offset += dom.GetMidDomainCollisionCount(2);

      StreamAndCollide(mOutletCollision, offset, dom.GetMidDomainCollisionCount(3));
      offset += dom.GetMidDomainCollisionCount(3);

      StreamAndCollide(mInletWallCollision, offset, dom.GetMidDomainCollisionCount(4));
      offset += dom.GetMidDomainCollisionCount(4);

      StreamAndCollide(mOutletWallCollision, offset, dom.GetMidDomainCollisionCount(5));

      timings[hemelb::reporting::Timers::lb_calc].Stop();
      timings[hemelb::reporting::Timers::lb].Stop();
    }

    template<class TRAITS>
    void LBM<TRAITS>::PostReceive()
    {
      timings[hemelb::reporting::Timers::lb].Start();

      // Copy the distribution functions received from the neighbouring
      // processors into the destination buffer "f_new".
      // This is done here, after receiving the sent distributions from neighbours.
      mLatDat->CopyReceived();

      auto& dom = mLatDat->GetDomain();
      // Do any cleanup steps necessary on boundary nodes
      site_t offset = dom.GetMidDomainSiteCount();

      timings[hemelb::reporting::Timers::lb_calc].Start();

      log::Logger::Log<log::Debug, log::OnePerCore>("LBM - PostReceive - StreamAndCollide");
      //TODO yup, this is horrible. If you read this, please improve the following code.
      PostStep(mMidFluidCollision, offset, dom.GetDomainEdgeCollisionCount(0));
      offset += dom.GetDomainEdgeCollisionCount(0);

      PostStep(mWallCollision, offset, dom.GetDomainEdgeCollisionCount(1));
      offset += dom.GetDomainEdgeCollisionCount(1);

      PostStep(mInletCollision, offset, dom.GetDomainEdgeCollisionCount(2));
      offset += dom.GetDomainEdgeCollisionCount(2);

      PostStep(mOutletCollision, offset, dom.GetDomainEdgeCollisionCount(3));
      offset += dom.GetDomainEdgeCollisionCount(3);

      PostStep(mInletWallCollision, offset, dom.GetDomainEdgeCollisionCount(4));
      offset += dom.GetDomainEdgeCollisionCount(4);

      PostStep(mOutletWallCollision, offset, dom.GetDomainEdgeCollisionCount(5));

      offset = 0;

      PostStep(mMidFluidCollision, offset, dom.GetMidDomainCollisionCount(0));
      offset += dom.GetMidDomainCollisionCount(0);

      PostStep(mWallCollision, offset, dom.GetMidDomainCollisionCount(1));
      offset += dom.GetMidDomainCollisionCount(1);

      PostStep(mInletCollision, offset, dom.GetMidDomainCollisionCount(2));
      offset += dom.GetMidDomainCollisionCount(2);

      PostStep(mOutletCollision, offset, dom.GetMidDomainCollisionCount(3));
      offset += dom.GetMidDomainCollisionCount(3);

      PostStep(mInletWallCollision, offset, dom.GetMidDomainCollisionCount(4));
      offset += dom.GetMidDomainCollisionCount(4);

      PostStep(mOutletWallCollision, offset, dom.GetMidDomainCollisionCount(5));

      timings[hemelb::reporting::Timers::lb_calc].Stop();
      timings[hemelb::reporting::Timers::lb].Stop();
    }

    template<class TRAITS>
    void LBM<TRAITS>::EndIteration()
    {
    }

    template<class TRAITS>
    LBM<TRAITS>::~LBM()
    {
      // Delete the collision and stream objects we've been using
      delete mMidFluidCollision;
      delete mWallCollision;
      delete mInletCollision;
      delete mOutletCollision;
      delete mInletWallCollision;
      delete mOutletWallCollision;
    }

  }
}

#endif /* HEMELB_LB_LB_HPP */
