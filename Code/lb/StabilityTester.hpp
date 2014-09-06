//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_LB_STABILITYTESTER_HPP
#define HEMELB_LB_STABILITYTESTER_HPP

namespace hemelb
{
  namespace lb
  {
    /**
     * Class to repeatedly assess the stability of the simulation, using non-blocking collective
     */
    template<class LatticeType>
    StabilityTester<LatticeType>::StabilityTester(const geometry::LatticeData * iLatDat,
                                                  net::Net* net, SimulationState* simState,
                                                  reporting::Timers& timings,
                                                 const hemelb::configuration::SimConfig::MonitoringConfig* testerConfig) :
        CollectiveAction(net->GetCommunicator(), timings[reporting::Timers::mpiWait]),
            mLatDat(iLatDat), mSimState(simState),
	    workTimer(timings[reporting::Timers::monitoring]), testerConfig(testerConfig)
    {
      Reset();
    }

    /**
     * Override the reset method in the base class, to reset the stability variables.
     */
    template<class LatticeType>
    void StabilityTester<LatticeType>::Reset()
    {
      localStability = UndefinedStability;
      globalStability = UndefinedStability;
      mSimState->SetStability(UndefinedStability);
    }

    /**
     * Compute the local stability/convergence state.
     */
    template<class LatticeType>
    void StabilityTester<LatticeType>::PreSend(void)
    {
      workTimer.Start();
      bool unconvergedSitePresent = false;
      bool checkConvThisTimeStep = testerConfig->doConvergenceCheck;
      localStability = Stable;

      for (site_t i = 0; i < mLatDat->GetLocalFluidSiteCount(); i++)
      {
        for (unsigned int l = 0; l < LatticeType::NUMVECTORS; l++)
        {
          distribn_t value = *mLatDat->GetFNew(i * LatticeType::NUMVECTORS + l);

          // Note that by testing for value > 0.0, we also catch stray NaNs.
          if (! (value > 0.0))
          {
            localStability = Unstable;
            break;
          }
        }
        ///@todo: If we refactor the previous loop out, we can get away with a single break statement
        if (localStability == Unstable)
          break;

        if (checkConvThisTimeStep)
        {
          distribn_t relativeDifference = ComputeRelativeDifference(mLatDat->GetFNew(i
                                                                        * LatticeType::NUMVECTORS),
                                                                    mLatDat->GetSite(i).GetFOld<
                                                                        LatticeType>());

          if (relativeDifference > testerConfig->convergenceRelativeTolerance)
          {
            // The simulation is stable but hasn't converged in the whole domain yet.
            unconvergedSitePresent = true;
            // Skip further tests for this timestep
            checkConvThisTimeStep = false;
          }
        }
      }

      if (localStability == Stable)
      {
        if (testerConfig->doConvergenceCheck && !unconvergedSitePresent)
        {
          localStability = StableAndConverged;
        }
      }
      workTimer.Stop();
    }

    /**
     * Initiate the collective.
     */
    template<class LatticeType>
    void StabilityTester<LatticeType>::Send(void)
    {
      // Begin collective.
      collectiveReq = collectiveComm.Iallreduce(localStability, MPI_MIN, globalStability);
    }

    /**
     * Computes the relative difference between the densities at the beginning and end of a
     * timestep, i.e. |(rho_new - rho_old) / (rho_old - rho_0)|.
     *
     * @param fNew Distribution function after stream and collide, i.e. solution of the current timestep.
     * @param fOld Distribution function at the end of the previous timestep.
     * @return relative difference between the densities computed from fNew and fOld.
     */
    template<class LatticeType>
    double StabilityTester<LatticeType>::ComputeRelativeDifference(const distribn_t* fNew,
                                                                   const distribn_t* fOld) const
    {
      distribn_t newDensity;
      distribn_t newMomentumX;
      distribn_t newMomentumY;
      distribn_t newMomentumZ;
      LatticeType::CalculateDensityAndMomentum(fNew,
					       newDensity,
					       newMomentumX,
					       newMomentumY,
					       newMomentumZ);
      
      distribn_t oldDensity;
      distribn_t oldMomentumX;
      distribn_t oldMomentumY;
      distribn_t oldMomentumZ;
      LatticeType::CalculateDensityAndMomentum(fOld,
					       oldDensity,
					       oldMomentumX,
					       oldMomentumY,
					       oldMomentumZ);
      
      distribn_t absoluteError;
      distribn_t referenceValue;
      
      switch (testerConfig->convergenceVariable)
      {
        case extraction::OutputField::Velocity:
	{
	  distribn_t diff_vel_x = newMomentumX / newDensity - oldMomentumX / oldDensity;
	  distribn_t diff_vel_y = newMomentumY / newDensity - oldMomentumY / oldDensity;
	  distribn_t diff_vel_z = newMomentumZ / newDensity - oldMomentumZ / oldDensity;
	  
	  absoluteError = sqrt(diff_vel_x * diff_vel_x + diff_vel_y * diff_vel_y
			       + diff_vel_z * diff_vel_z);
	  referenceValue = testerConfig->convergenceReferenceValue;
	  break;
	}
        default:
	  // Never reached
	  throw Exception()
	    << "Convergence check based on requested variable currently not available";
      }
      
      return absoluteError / referenceValue;
    }

    /**
     * Apply the stability value sent by the root node to the simulation logic.
     */
    template<class LatticeType>
    void StabilityTester<LatticeType>::PostReceive()
    {
      mSimState->SetStability((Stability) globalStability);
    }
  }
}

#endif /* HEMELB_LB_STABILITYTESTER_HPP */
