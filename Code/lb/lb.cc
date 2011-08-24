// In this file, the functions useful to calculate the equilibrium distribution
// function, momentums, the effective von Mises stress and the boundary conditions
// are reported

#include <math.h>
#include <limits>

#include "lb/lb.h"
#include "util/utilityFunctions.h"
#include "vis/RayTracer.h"
#include "lb/rheology_models/CassonRheologyModel.h"
#include "lb/rheology_models/TruncatedPowerLawRheologyModel.h"

namespace hemelb
{
  namespace lb
  {
    distribn_t LBM::ConvertPressureToLatticeUnits(double pressure) const
    {
      double temp = (PULSATILE_PERIOD_s / ((double) mState->GetTimeStepsPerCycle() * mLatDat->GetVoxelSize()));
      return Cs2
          + (pressure - REFERENCE_PRESSURE_mmHg) * mmHg_TO_PASCAL * temp * temp
              / BLOOD_DENSITY_Kg_per_m3;
    }

    double LBM::ConvertPressureToPhysicalUnits(distribn_t pressure) const
    {
      double temp = ( ((double) mState->GetTimeStepsPerCycle() * mLatDat->GetVoxelSize()) / PULSATILE_PERIOD_s);
      return REFERENCE_PRESSURE_mmHg
          + ( (pressure / Cs2 - 1.0) * Cs2) * BLOOD_DENSITY_Kg_per_m3 * temp * temp / mmHg_TO_PASCAL;
    }

    distribn_t LBM::ConvertPressureGradToLatticeUnits(double pressure_grad) const
    {
      double temp = (PULSATILE_PERIOD_s / ((double) mState->GetTimeStepsPerCycle() * mLatDat->GetVoxelSize()));
      return pressure_grad * mmHg_TO_PASCAL * temp * temp / BLOOD_DENSITY_Kg_per_m3;
    }

    double LBM::ConvertPressureGradToPhysicalUnits(distribn_t pressure_grad) const
    {
      double temp = ( ((double) mState->GetTimeStepsPerCycle() * mLatDat->GetVoxelSize()) / PULSATILE_PERIOD_s);
      return pressure_grad * BLOOD_DENSITY_Kg_per_m3 * temp * temp / mmHg_TO_PASCAL;
    }

    distribn_t LBM::ConvertVelocityToLatticeUnits(double velocity) const
    {
      return velocity * ( ( (mParams.Tau - 0.5) / 3.0) * mLatDat->GetVoxelSize())
          / (BLOOD_VISCOSITY_Pa_s / BLOOD_DENSITY_Kg_per_m3);
    }

    double LBM::ConvertVelocityToPhysicalUnits(distribn_t velocity) const
    {
      // convert velocity from lattice units to physical units (m/s)
      return velocity * (BLOOD_VISCOSITY_Pa_s / BLOOD_DENSITY_Kg_per_m3)
          / ( ( (mParams.Tau - 0.5) / 3.0) * mLatDat->GetVoxelSize());
    }

    distribn_t LBM::ConvertStressToLatticeUnits(double stress) const
    {
      return stress * (BLOOD_DENSITY_Kg_per_m3 / (BLOOD_VISCOSITY_Pa_s * BLOOD_VISCOSITY_Pa_s))
          * ( ( (mParams.Tau - 0.5) / 3.0) * mLatDat->GetVoxelSize())
          * ( ( (mParams.Tau - 0.5) / 3.0) * mLatDat->GetVoxelSize());
    }

    double LBM::ConvertStressToPhysicalUnits(distribn_t stress) const
    {
      // convert stress from lattice units to physical units (Pa)
      return stress * BLOOD_VISCOSITY_Pa_s * BLOOD_VISCOSITY_Pa_s
          / (BLOOD_DENSITY_Kg_per_m3 * ( ( (mParams.Tau - 0.5) / 3.0) * mLatDat->GetVoxelSize())
              * ( ( (mParams.Tau - 0.5) / 3.0) * mLatDat->GetVoxelSize()));
    }

    void LBM::RecalculateTauViscosityOmega()
    {
      mParams.Tau = 0.5
          + (PULSATILE_PERIOD_s * BLOOD_VISCOSITY_Pa_s / BLOOD_DENSITY_Kg_per_m3)
              / (Cs2 * ((double) mState->GetTimeStepsPerCycle() * mLatDat->GetVoxelSize() * mLatDat->GetVoxelSize()));

      mParams.Omega = -1.0 / mParams.Tau;
      mParams.StressParameter = (1.0 - 1.0 / (2.0 * mParams.Tau)) / sqrt(2.0);
      mParams.Beta = -1.0 / (2.0 * mParams.Tau);
    }

    // Calculate the BCs for each boundary site type and the
    // non-equilibrium distribution functions.
    void LBM::CalculateBC(distribn_t f[],
                          hemelb::geometry::LatticeData::SiteType iSiteType,
                          unsigned int iBoundaryId,
                          distribn_t *density,
                          distribn_t *vx,
                          distribn_t *vy,
                          distribn_t *vz,
                          distribn_t f_neq[])
    {
      distribn_t dummy_density;

      for (unsigned int l = 0; l < D3Q15::NUMVECTORS; l++)
      {
        f_neq[l] = f[l];
      }

      if (iSiteType == hemelb::geometry::LatticeData::FLUID_TYPE)
      {
        D3Q15::CalculateDensityAndVelocity(f, *density, *vx, *vy, *vz);
      }
      else
      {
        if (iSiteType == hemelb::geometry::LatticeData::INLET_TYPE)
        {
          *density = inlet_density[iBoundaryId];
        }
        else
        {
          *density = outlet_density[iBoundaryId];
        }

        D3Q15::CalculateDensityAndVelocity(f, dummy_density, *vx, *vy, *vz);
        D3Q15::CalculateFeq(*density, *vx, *vy, *vz, f);

      }
      for (unsigned int l = 0; l < D3Q15::NUMVECTORS; l++)
      {
        f_neq[l] -= f[l];
      }

    }

    void LBM::UpdateBoundaryDensities(unsigned long time_step)
    {
      double w = 2.0 * PI / (double) mState->GetTimeStepsPerCycle();

      for (int i = 0; i < inlets; i++)
      {
        inlet_density[i] = inlet_density_avg[i]
            + inlet_density_amp[i] * cos(w * (double) time_step + inlet_density_phs[i]);
      }
      for (int i = 0; i < outlets; i++)
      {
        outlet_density[i] = outlet_density_avg[i]
            + outlet_density_amp[i] * cos(w * (double) time_step + outlet_density_phs[i]);
      }
    }

    hemelb::lb::LbmParameters *LBM::GetLbmParams()
    {
      return &mParams;
    }

    LBM::LBM(SimConfig *iSimulationConfig,
             net::Net* net,
             geometry::LatticeData* latDat,
             SimulationState* simState) :
        mSimConfig(iSimulationConfig), mNet(net), mLatDat(latDat), mState(simState)
    {
      // voxel_size = iSimulationConfig->VoxelSize;

      ReadParameters();

      //typedef hemelb::lb::collisions::implementations::LBGK CO;

      //typedef hemelb::lb::collisions::implementations::LBGKNN<hemelb::lb::rheology_models::CassonRheologyModel> CO;
      typedef hemelb::lb::collisions::implementations::LBGKNN<hemelb::lb::rheology_models::TruncatedPowerLawRheologyModel> CO;

      if (typeid(CO) == typeid(hemelb::lb::collisions::implementations::ELBM))
      {
        hemelb::lb::collisions::implementations::ELBM::createAlphaArray(mLatDat->GetLocalFluidSiteCount());
        hemelb::lb::collisions::implementations::ELBM::setTau(&mParams.Tau);
      }

      if (typeid(CO) == typeid(hemelb::lb::collisions::implementations::LBGKNN<hemelb::lb::rheology_models::CassonRheologyModel>))
      {
        hemelb::lb::collisions::implementations::LBGKNN<hemelb::lb::rheology_models::CassonRheologyModel>::createTauArray(mLatDat->GetLocalFluidSiteCount(),
                                                                                                                          mParams.Tau);
        hemelb::lb::collisions::implementations::LBGKNN<hemelb::lb::rheology_models::CassonRheologyModel>::setStateObjects(mLatDat, mState);
      }

      if (typeid(CO) == typeid(hemelb::lb::collisions::implementations::LBGKNN<hemelb::lb::rheology_models::TruncatedPowerLawRheologyModel>))
      {
        hemelb::lb::collisions::implementations::LBGKNN<hemelb::lb::rheology_models::TruncatedPowerLawRheologyModel>::createTauArray(mLatDat->GetLocalFluidSiteCount(),
                                                                                                                                     mParams.Tau);
        hemelb::lb::collisions::implementations::LBGKNN<hemelb::lb::rheology_models::TruncatedPowerLawRheologyModel>::setStateObjects(mLatDat, mState);
      }


      InitCollisions<hemelb::lb::streamers::implementations::SimpleCollideAndStream<CO>,
          hemelb::lb::streamers::implementations::ZeroVelocityEquilibrium<CO>,
          hemelb::lb::streamers::implementations::NonZeroVelocityBoundaryDensity<CO>,
          hemelb::lb::streamers::implementations::ZeroVelocityBoundaryDensity<CO> >();
    }

    void LBM::CalculateMouseFlowField(float densityIn,
                                      float stressIn,
                                      distribn_t &mouse_pressure,
                                      distribn_t &mouse_stress,
                                      double density_threshold_min,
                                      double density_threshold_minmax_inv,
                                      double stress_threshold_max_inv)
    {
      double density = density_threshold_min + densityIn / density_threshold_minmax_inv;
      double stress = stressIn / stress_threshold_max_inv;

      mouse_pressure = ConvertPressureToPhysicalUnits(density * Cs2);
      mouse_stress = ConvertStressToPhysicalUnits(stress);
    }

    template<typename tMidFluidCollision, typename tWallCollision, typename tInletOutletCollision,
        typename tInletOutletWallCollision>
    void LBM::InitCollisions()
    {
      mStreamAndCollide = new hemelb::lb::collisions::StreamAndCollide<tMidFluidCollision,
          tWallCollision, tInletOutletCollision, tInletOutletWallCollision>();
      mPostStep = new hemelb::lb::collisions::PostStep<tMidFluidCollision, tWallCollision,
          tInletOutletCollision, tInletOutletWallCollision>();

      // TODO Note that the convergence checking is not yet implemented in the
      // new boundary condition hierarchy system.
      // It'd be nice to do this with something like
      // MidFluidCollision = new ConvergenceCheckingWrapper(new WhateverMidFluidCollision());

      mMidFluidCollision = new hemelb::lb::streamers::MidFluidCollision();
      mWallCollision = new hemelb::lb::streamers::WallCollision();
      mInletCollision = new hemelb::lb::streamers::InletOutletCollision(inlet_density);
      mOutletCollision = new hemelb::lb::streamers::InletOutletCollision(outlet_density);
      mInletWallCollision = new hemelb::lb::streamers::InletOutletWallCollision(inlet_density);
      mOutletWallCollision = new hemelb::lb::streamers::InletOutletWallCollision(outlet_density);
    }

    void LBM::Initialise(site_t* iFTranslator, vis::Control* iControl)
    {
      receivedFTranslator = iFTranslator;

      SetInitialConditions();

      mVisControl = iControl;
    }

    void LBM::SetInitialConditions()
    {
      timeSpent = 0.0;

      distribn_t *f_old_p, *f_new_p, f_eq[D3Q15::NUMVECTORS];
      distribn_t density;

      density = 0.;

      for (int i = 0; i < outlets; i++)
      {
        density += outlet_density_avg[i] - outlet_density_amp[i];
      }
      density /= outlets;

      for (site_t i = 0; i < mLatDat->GetLocalFluidSiteCount(); i++)
      {
        D3Q15::CalculateFeq(density, 0.0, 0.0, 0.0, f_eq);

        f_old_p = mLatDat->GetFOld(i * D3Q15::NUMVECTORS);
        f_new_p = mLatDat->GetFNew(i * D3Q15::NUMVECTORS);

        for (unsigned int l = 0; l < D3Q15::NUMVECTORS; l++)
        {
          f_new_p[l] = f_old_p[l] = f_eq[l];
        }
      }
    }

    // TODO HACK
    hemelb::lb::streamers::Collision* LBM::GetCollision(int i)
    {
      switch (i)
      {
        case 0:
          return mMidFluidCollision;
        case 1:
          return mWallCollision;
        case 2:
          return mInletCollision;
        case 3:
          return mOutletCollision;
        case 4:
          return mInletWallCollision;
        case 5:
          return mOutletWallCollision;
      }
      return NULL;
    }

    void LBM::RequestComms()
    {
      timeSpent -= util::myClock();

      topology::NetworkTopology* netTop = topology::NetworkTopology::Instance();

      for (std::vector<hemelb::topology::NeighbouringProcessor>::const_iterator it =
          netTop->NeighbouringProcs.begin(); it != netTop->NeighbouringProcs.end(); it++)
      {
        // Request the receive into the appropriate bit of FOld.
        mNet->RequestReceive<distribn_t>(mLatDat->GetFOld( (*it).FirstSharedF),
                                         (int) (*it).SharedFCount,
                                         (*it).Rank);

        // Request the send from the right bit of FNew.
        mNet->RequestSend<distribn_t>(mLatDat->GetFNew( (*it).FirstSharedF),
                                      (int) (*it).SharedFCount,
                                      (*it).Rank);
      }

      timeSpent += util::myClock();
    }

    void LBM::PreSend()
    {
      timeSpent -= util::myClock();

      site_t offset = mLatDat->GetInnerSiteCount();

      for (unsigned int collision_type = 0; collision_type < COLLISION_TYPES; collision_type++)
      {
        GetCollision(collision_type)->AcceptCollisionVisitor(mStreamAndCollide,
                                                             mVisControl->IsRendering(),
                                                             offset,
                                                             mLatDat->GetInterCollisionCount(collision_type),
                                                             &mParams,
                                                             mLatDat,
                                                             mVisControl);
        offset += mLatDat->GetInterCollisionCount(collision_type);
      }

      timeSpent += util::myClock();
    }

    void LBM::PreReceive()
    {
      timeSpent -= util::myClock();

      site_t offset = 0;

      for (unsigned int collision_type = 0; collision_type < COLLISION_TYPES; collision_type++)
      {
        GetCollision(collision_type)->AcceptCollisionVisitor(mStreamAndCollide,
                                                             mVisControl->IsRendering(),
                                                             offset,
                                                             mLatDat->GetInnerCollisionCount(collision_type),
                                                             &mParams,
                                                             mLatDat,
                                                             mVisControl);
        offset += mLatDat->GetInnerCollisionCount(collision_type);
      }

      timeSpent += util::myClock();
    }

    void LBM::PostReceive()
    {
      timeSpent -= util::myClock();

      // Copy the distribution functions received from the neighbouring
      // processors into the destination buffer "f_new".
      topology::NetworkTopology* netTop = topology::NetworkTopology::Instance();

      for (site_t i = 0; i < netTop->TotalSharedFs; i++)
      {
        *mLatDat->GetFNew(receivedFTranslator[i]) =
            *mLatDat->GetFOld(netTop->NeighbouringProcs[0].FirstSharedF + i);
      }

      // Do any cleanup steps necessary on boundary nodes
      size_t offset = 0;

      for (unsigned int collision_type = 0; collision_type < COLLISION_TYPES; collision_type++)
      {
        GetCollision(collision_type)->AcceptCollisionVisitor(mPostStep,
                                                             mVisControl->IsRendering(),
                                                             offset,
                                                             mLatDat->GetInnerCollisionCount(collision_type),
                                                             &mParams,
                                                             mLatDat,
                                                             mVisControl);
        offset += mLatDat->GetInnerCollisionCount(collision_type);
      }

      for (unsigned int collision_type = 0; collision_type < COLLISION_TYPES; collision_type++)
      {
        GetCollision(collision_type)->AcceptCollisionVisitor(mPostStep,
                                                             mVisControl->IsRendering(),
                                                             offset,
                                                             mLatDat->GetInterCollisionCount(collision_type),
                                                             &mParams,
                                                             mLatDat,
                                                             mVisControl);
        offset += mLatDat->GetInterCollisionCount(collision_type);
      }

      timeSpent += util::myClock();
    }

    void LBM::EndIteration()
    {
      timeSpent -= util::myClock();

      // Swap f_old and f_new ready for the next timestep.
      mLatDat->SwapOldAndNew();

      timeSpent += util::myClock();
    }

    // Update peak and average inlet velocities local to the current subdomain.
    void LBM::UpdateInletVelocities(unsigned long time_step)
    {
      distribn_t density;
      distribn_t vx, vy, vz;
      distribn_t velocity;

      int inlet_id;

      site_t offset = mLatDat->GetInnerCollisionCount(0) + mLatDat->GetInnerCollisionCount(1);

      for (site_t i = offset; i < offset + mLatDat->GetInnerCollisionCount(2); i++)
      {
        D3Q15::CalculateDensityAndVelocity(mLatDat->GetFOld(i * D3Q15::NUMVECTORS),
                                           density,
                                           vx,
                                           vy,
                                           vz);

        inlet_id = mLatDat->GetBoundaryId(i);

        vx *= inlet_normal[3 * inlet_id + 0];
        vy *= inlet_normal[3 * inlet_id + 1];
        vz *= inlet_normal[3 * inlet_id + 2];

        velocity = vx * vx + vy * vy + vz * vz;

        if (velocity > 0.)
        {
          velocity = sqrt(velocity) / density;
        }
        else
        {
          velocity = -sqrt(velocity) / density;
        }
      }

      offset = mLatDat->GetInnerSiteCount() + mLatDat->GetInterCollisionCount(0)
          + mLatDat->GetInterCollisionCount(1);

      for (site_t i = offset; i < offset + mLatDat->GetInterCollisionCount(2); i++)
      {
        D3Q15::CalculateDensityAndVelocity(mLatDat->GetFOld(i * D3Q15::NUMVECTORS),
                                           density,
                                           vx,
                                           vy,
                                           vz);

        inlet_id = mLatDat->GetBoundaryId(i);

        vx *= inlet_normal[3 * inlet_id + 0];
        vy *= inlet_normal[3 * inlet_id + 1];
        vz *= inlet_normal[3 * inlet_id + 2];

        velocity = vx * vx + vy * vy + vz * vz;

        if (velocity > 0.)
        {
          velocity = sqrt(velocity) / density;
        }
        else
        {
          velocity = -sqrt(velocity) / density;
        }
      }
    }

    /**
     * Return the amount of time spent doing lattice-Boltzmann
     * @return
     */
    double LBM::GetTimeSpent() const
    {
      return timeSpent;
    }

    // In the case of instability, this function restart the simulation
    // with twice as many time steps per period and update the parameters
    // that depends on this change.
    void LBM::Reset()
    {
      int i;

      for (i = 0; i < inlets; i++)
      {
        inlet_density_avg[i] = ConvertPressureToPhysicalUnits(inlet_density_avg[i] * Cs2);
        inlet_density_amp[i] = ConvertPressureGradToPhysicalUnits(inlet_density_amp[i] * Cs2);
      }
      for (i = 0; i < outlets; i++)
      {
        outlet_density_avg[i] = ConvertPressureToPhysicalUnits(outlet_density_avg[i] * Cs2);
        outlet_density_amp[i] = ConvertPressureGradToPhysicalUnits(outlet_density_amp[i] * Cs2);
      }

      mState->DoubleTimeResolution();

      for (i = 0; i < inlets; i++)
      {
        inlet_density_avg[i] = ConvertPressureToLatticeUnits(inlet_density_avg[i]) / Cs2;
        inlet_density_amp[i] = ConvertPressureGradToLatticeUnits(inlet_density_amp[i]) / Cs2;
      }
      for (i = 0; i < outlets; i++)
      {
        outlet_density_avg[i] = ConvertPressureToLatticeUnits(outlet_density_avg[i]) / Cs2;
        outlet_density_amp[i] = ConvertPressureGradToLatticeUnits(outlet_density_amp[i]) / Cs2;
      }

      RecalculateTauViscosityOmega();

      SetInitialConditions();
    }

    LBM::~LBM()
    {
      // Delete the translator between received location and location in f_new.
      delete[] receivedFTranslator;

      // Delete arrays allocated for the inlets
      delete[] inlet_density;
      delete[] inlet_density_avg;
      delete[] inlet_density_amp;
      delete[] inlet_density_phs;

      // Delete arrays allocated for the outlets
      delete[] outlet_density;
      delete[] outlet_density_avg;
      delete[] outlet_density_amp;
      delete[] outlet_density_phs;

      // Delete visitors
      delete mStreamAndCollide;
      delete mPostStep;

      // Delete the collision and stream objects we've been using
      delete mMidFluidCollision;
      delete mWallCollision;
      delete mInletCollision;
      delete mOutletCollision;
      delete mInletWallCollision;
      delete mOutletWallCollision;

      // Delete various other arrays used
      delete[] inlet_normal;
    }
  }
}
