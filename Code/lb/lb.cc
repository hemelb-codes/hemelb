// In this file, the functions useful to calculate the equilibrium distribution
// function, momentums, the effective von Mises stress and the boundary conditions
// are reported

#include <math.h>
#include <limits>

#include "lb/lb.h"
#include "util/utilityFunctions.h"
#include "vis/RayTracer.h"

namespace hemelb
{
  namespace lb
  {
    double LBM::ConvertPressureToLatticeUnits(double pressure) const
    {
      return Cs2 + (pressure - REFERENCE_PRESSURE) * mmHg_TO_PASCAL * (PULSATILE_PERIOD / (period
          * voxel_size)) * (PULSATILE_PERIOD / (period * voxel_size)) / BLOOD_DENSITY;
    }

    double LBM::ConvertPressureToPhysicalUnits(double pressure) const
    {
      return REFERENCE_PRESSURE + ( (pressure / Cs2 - 1.0) * Cs2) * BLOOD_DENSITY * ( (period
          * voxel_size) / PULSATILE_PERIOD) * ( (period * voxel_size) / PULSATILE_PERIOD)
          / mmHg_TO_PASCAL;
    }

    double LBM::ConvertPressureGradToLatticeUnits(double pressure_grad) const
    {
      return pressure_grad * mmHg_TO_PASCAL * (PULSATILE_PERIOD / (period * voxel_size))
          * (PULSATILE_PERIOD / (period * voxel_size)) / BLOOD_DENSITY;
    }

    double LBM::ConvertPressureGradToPhysicalUnits(double pressure_grad) const
    {
      return pressure_grad * BLOOD_DENSITY * ( (period * voxel_size) / PULSATILE_PERIOD)
          * ( (period * voxel_size) / PULSATILE_PERIOD) / mmHg_TO_PASCAL;
    }

    double LBM::ConvertVelocityToLatticeUnits(double velocity) const
    {
      return velocity * ( ( (mParams.Tau - 0.5) / 3.0) * voxel_size) / (BLOOD_VISCOSITY
          / BLOOD_DENSITY);
    }

    double LBM::ConvertVelocityToPhysicalUnits(double velocity) const
    {
      // convert velocity from lattice units to physical units (m/s)
      return velocity * (BLOOD_VISCOSITY / BLOOD_DENSITY) / ( ( (mParams.Tau - 0.5) / 3.0)
          * voxel_size);
    }

    double LBM::ConvertStressToLatticeUnits(double stress) const
    {
      return stress * (BLOOD_DENSITY / (BLOOD_VISCOSITY * BLOOD_VISCOSITY)) * ( ( (mParams.Tau
          - 0.5) / 3.0) * voxel_size) * ( ( (mParams.Tau - 0.5) / 3.0) * voxel_size);
    }

    double LBM::ConvertStressToPhysicalUnits(double stress) const
    {
      // convert stress from lattice units to physical units (Pa)
      return stress * BLOOD_VISCOSITY * BLOOD_VISCOSITY / (BLOOD_DENSITY * ( ( (mParams.Tau - 0.5)
          / 3.0) * voxel_size) * ( ( (mParams.Tau - 0.5) / 3.0) * voxel_size));
    }

    void LBM::RecalculateTauViscosityOmega()
    {
      mParams.Tau = 0.5 + (PULSATILE_PERIOD * BLOOD_VISCOSITY / BLOOD_DENSITY) / (Cs2 * period
          * voxel_size * voxel_size);

      mParams.Omega = -1.0 / mParams.Tau;
      mParams.StressParameter = (1.0 - 1.0 / (2.0 * mParams.Tau)) / sqrt(2.0);
    }

    // Calculate the BCs for each boundary site type and the
    // non-equilibrium distribution functions.
    void LBM::CalculateBC(double f[],
                          hemelb::geometry::LatticeData::SiteType iSiteType,
                          unsigned int iBoundaryId,
                          double *density,
                          double *vx,
                          double *vy,
                          double *vz,
                          double f_neq[])
    {
      double dummy_density;

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

    void LBM::UpdateBoundaryDensities(int cycle_id, int time_step)
    {
      double w = 2.0 * PI / period;

      for (int i = 0; i < inlets; i++)
      {
        inlet_density[i] = inlet_density_avg[i] + inlet_density_amp[i] * cos(w * (double) time_step
            + inlet_density_phs[i]);
      }
      for (int i = 0; i < outlets; i++)
      {
        outlet_density[i] = outlet_density_avg[i] + outlet_density_amp[i] * cos(w
            * (double) time_step + outlet_density_phs[i]);
      }
    }

    hemelb::lb::LbmParameters *LBM::GetLbmParams()
    {
      return &mParams;
    }

    LBM::LBM(SimConfig *iSimulationConfig,
             net::Net* net,
             geometry::LatticeData* latDat,
             SimulationState* simState,
             const topology::NetworkTopology * iNetTop) :
      mSimConfig(iSimulationConfig), mNet(net), mLatDat(latDat), mState(simState),
          mNetTopology(iNetTop)
    {
      period = iSimulationConfig->StepsPerCycle;
      voxel_size = iSimulationConfig->VoxelSize;

      ReadParameters();

      InitCollisions();
    }

    void LBM::CalculateMouseFlowField(float densityIn,
                                      float stressIn,
                                      double &mouse_pressure,
                                      double &mouse_stress,
                                      double density_threshold_min,
                                      double density_threshold_minmax_inv,
                                      double stress_threshold_max_inv)
    {
      double density = density_threshold_min + densityIn / density_threshold_minmax_inv;
      double stress = stressIn / stress_threshold_max_inv;

      mouse_pressure = ConvertPressureToPhysicalUnits(density * Cs2);
      mouse_stress = ConvertStressToPhysicalUnits(stress);
    }

    void LBM::InitCollisions()
    {
      // TODO Note that the convergence checking is not yet implemented in the
      // new boundary condition hierarchy system.
      // It'd be nice to do this with something like
      // MidFluidCollision = new ConvergenceCheckingWrapper(new WhateverMidFluidCollision());

      mMidFluidCollision = new hemelb::lb::collisions::ImplSimpleCollideAndStream();
      mWallCollision = new hemelb::lb::collisions::ImplZeroVelocityEquilibrium();
      mInletCollision
          = new hemelb::lb::collisions::ImplNonZeroVelocityBoundaryDensity(inlet_density);
      mOutletCollision
          = new hemelb::lb::collisions::ImplNonZeroVelocityBoundaryDensity(outlet_density);
      mInletWallCollision
          = new hemelb::lb::collisions::ImplZeroVelocityBoundaryDensity(inlet_density);
      mOutletWallCollision
          = new hemelb::lb::collisions::ImplZeroVelocityBoundaryDensity(outlet_density);
    }

    void LBM::Initialise(int* iFTranslator, vis::Control* iControl)
    {
      receivedFTranslator = iFTranslator;

      SetInitialConditions();

      mVisControl = iControl;
    }

    void LBM::SetInitialConditions()
    {
      double *f_old_p, *f_new_p, f_eq[D3Q15::NUMVECTORS];
      double density;

      density = 0.;

      for (int i = 0; i < outlets; i++)
      {
        density += outlet_density_avg[i] - outlet_density_amp[i];
      }
      density /= outlets;

      for (unsigned int i = 0; i < mLatDat->GetLocalFluidSiteCount(); i++)
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
    hemelb::lb::collisions::Collision* LBM::GetCollision(int i)
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
      for (std::vector<hemelb::topology::NeighbouringProcessor*>::const_iterator it =
          mNetTopology->NeighbouringProcs.begin(); it != mNetTopology->NeighbouringProcs.end(); it++)
      {
        // Request the receive into the appropriate bit of FOld.
        mNet->RequestReceive<double> (mLatDat->GetFOld( (*it)->FirstSharedF),
                                       (*it)->SharedFCount,
                                       (*it)->Rank);

        // Request the send from the right bit of
        mNet->RequestSend<double> (mLatDat->GetFNew( (*it)->FirstSharedF),
                                    (*it)->SharedFCount,
                                    (*it)->Rank);
      }
    }

    void LBM::PreSend()
    {
      int offset = mLatDat->GetInnerSiteCount();

      for (int collision_type = 0; collision_type < COLLISION_TYPES; collision_type++)
      {
        GetCollision(collision_type)->DoCollisions(mState->DoRendering,
                                                   offset,
                                                   mLatDat->GetInterCollisionCount(collision_type),
                                                   mParams,
                                                   mMinsAndMaxes,
                                                   *mLatDat,
                                                   mVisControl);
        offset += mLatDat->GetInterCollisionCount(collision_type);
      }
    }

    void LBM::PreReceive()
    {
      int offset = 0;

      for (int collision_type = 0; collision_type < COLLISION_TYPES; collision_type++)
      {
        GetCollision(collision_type)->DoCollisions(mState->DoRendering,
                                                   offset,
                                                   mLatDat->GetInnerCollisionCount(collision_type),
                                                   mParams,
                                                   mMinsAndMaxes,
                                                   *mLatDat,
                                                   mVisControl);
        offset += mLatDat->GetInnerCollisionCount(collision_type);
      }
    }

    void LBM::PostReceive()
    {
      // Copy the distribution functions received from the neighbouring
      // processors into the destination buffer "f_new".
      for (int i = 0; i < mNetTopology->TotalSharedFs; i++)
      {
        *mLatDat->GetFNew(receivedFTranslator[i])
            = *mLatDat->GetFOld(mNetTopology->NeighbouringProcs[0]->FirstSharedF + i);
      }

      // Do any cleanup steps necessary on boundary nodes
      int offset = 0;

      for (int collision_type = 0; collision_type < COLLISION_TYPES; collision_type++)
      {
        GetCollision(collision_type)->PostStep(mState->DoRendering,
                                               offset,
                                               mLatDat->GetInnerCollisionCount(collision_type),
                                               mParams,
                                               mMinsAndMaxes,
                                               *mLatDat,
                                               mVisControl);
        offset += mLatDat->GetInnerCollisionCount(collision_type);
      }

      for (int collision_type = 0; collision_type < COLLISION_TYPES; collision_type++)
      {
        GetCollision(collision_type)->PostStep(mState->DoRendering,
                                               offset,
                                               mLatDat->GetInterCollisionCount(collision_type),
                                               mParams,
                                               mMinsAndMaxes,
                                               *mLatDat,
                                               mVisControl);
        offset += mLatDat->GetInterCollisionCount(collision_type);
      }
    }

    void LBM::EndIteration()
    {
      // Swap f_old and f_new ready for the next timestep.
      mLatDat->SwapOldAndNew();
    }

    void LBM::CalculateFlowFieldValues()
    {
      double *local_data;
      double *global_data;

      int i;

      int lMaxInlets = hemelb::util::NumericalFunctions::max<int>(6 + inlets, 2 * inlets);

      local_data = new double[lMaxInlets];
      global_data = new double[lMaxInlets];

      local_data[0] = mMinsAndMaxes.MinDensity;
      local_data[1] = mMinsAndMaxes.MinVelocity;
      local_data[2] = mMinsAndMaxes.MinStress;
      local_data[3] = mMinsAndMaxes.MaxDensity;
      local_data[4] = mMinsAndMaxes.MaxVelocity;
      local_data[5] = mMinsAndMaxes.MaxStress;

      memcpy(&local_data[6], peak_inlet_velocity, sizeof(double) * inlets);

      MPI_Reduce(&local_data[0], &global_data[0], 3, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
      MPI_Reduce(&local_data[3],
                 &global_data[3],
                 3 + inlets,
                 MPI_DOUBLE,
                 MPI_MAX,
                 0,
                 MPI_COMM_WORLD);

      mMinsAndMaxes.MinDensity = global_data[0];
      mMinsAndMaxes.MinVelocity = global_data[1];
      mMinsAndMaxes.MinStress = global_data[2];
      mMinsAndMaxes.MaxDensity = global_data[3];
      mMinsAndMaxes.MaxVelocity = global_data[4];
      mMinsAndMaxes.MaxStress = global_data[5];

      memcpy(peak_inlet_velocity, &global_data[6], sizeof(double) * inlets);

      for (i = 0; i < inlets; i++)
      {
        local_data[i] = average_inlet_velocity[i];
        local_data[inlets + i] = inlet_count[i];
      }
      MPI_Reduce(local_data, global_data, 2 * inlets, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

      for (i = 0; i < inlets; i++)
      {
        average_inlet_velocity[i] = global_data[i];
        inlet_count[i] = global_data[inlets + i];
      }

      delete[] global_data;
      delete[] local_data;

      for (i = 0; i < inlets; i++)
      {
        average_inlet_velocity[i] /= inlet_count[i];
        average_inlet_velocity[i] = ConvertVelocityToPhysicalUnits(average_inlet_velocity[i]);
        peak_inlet_velocity[i] = ConvertVelocityToPhysicalUnits(peak_inlet_velocity[i]);
      }
    }

    // Update peak and average inlet velocities local to the current subdomain.
    void LBM::UpdateInletVelocities(int time_step)
    {
      double density;
      double vx, vy, vz;
      double velocity;

      int inlet_id;
      int c1, c2;

      if (time_step == 1)
      {
        for (int i = 0; i < inlets; i++)
        {
          peak_inlet_velocity[i] = -BIG_NUMBER;
          average_inlet_velocity[i] = 0.;
          inlet_count[i] = 0;
        }
      }

      unsigned int offset = mLatDat->GetInnerCollisionCount(0) + mLatDat->GetInnerCollisionCount(1);

      for (unsigned int i = offset; i < offset + mLatDat->GetInnerCollisionCount(2); i++)
      {
        D3Q15::CalculateDensityAndVelocity(mLatDat->GetFOld(i * D3Q15::NUMVECTORS),
                                           density,
                                           vx,
                                           vy,
                                           vz);

        inlet_id = mLatDat->GetBoundaryId(i);

        if (is_inlet_normal_available)
        {
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
        else
        {
          velocity = sqrt(vx * vx + vy * vy + vz * vz) / density;
        }
        peak_inlet_velocity[inlet_id] = fmax(peak_inlet_velocity[inlet_id], velocity);
        average_inlet_velocity[inlet_id] += velocity;
        ++inlet_count[inlet_id];
      }

      offset = mLatDat->GetInnerSiteCount() + mLatDat->GetInterCollisionCount(0)
          + mLatDat->GetInterCollisionCount(1);

      for (unsigned int i = offset; i < offset + mLatDat->GetInterCollisionCount(2); i++)
      {
        D3Q15::CalculateDensityAndVelocity(mLatDat->GetFOld(i * D3Q15::NUMVECTORS),
                                           density,
                                           vx,
                                           vy,
                                           vz);

        inlet_id = mLatDat->GetBoundaryId(i);

        if (is_inlet_normal_available)
        {
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
        else
        {
          velocity = sqrt(vx * vx + vy * vy + vz * vz) / density;
        }
        peak_inlet_velocity[inlet_id] = fmax(peak_inlet_velocity[inlet_id], velocity);
        average_inlet_velocity[inlet_id] += velocity;
        ++inlet_count[inlet_id];
      }
    }

    double LBM::GetAverageInletVelocity(int iInletNumber)
    {
      return average_inlet_velocity[iInletNumber];
    }
    double LBM::GetPeakInletVelocity(int iInletNumber)
    {
      return peak_inlet_velocity[iInletNumber];
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
      period *= 2;

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

    double LBM::GetMinPhysicalPressure()
    {
      return ConvertPressureToPhysicalUnits(mMinsAndMaxes.MinDensity * Cs2);
    }
    double LBM::GetMaxPhysicalPressure()
    {
      return ConvertPressureToPhysicalUnits(mMinsAndMaxes.MaxDensity * Cs2);
    }
    double LBM::GetMinPhysicalVelocity()
    {
      return ConvertVelocityToPhysicalUnits(mMinsAndMaxes.MinVelocity);
    }
    double LBM::GetMaxPhysicalVelocity()
    {
      return ConvertVelocityToPhysicalUnits(mMinsAndMaxes.MaxVelocity);
    }
    double LBM::GetMinPhysicalStress()
    {
      return ConvertStressToPhysicalUnits(mMinsAndMaxes.MinStress);
    }
    double LBM::GetMaxPhysicalStress()
    {
      return ConvertStressToPhysicalUnits(mMinsAndMaxes.MaxStress);
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

      // Delete the collision and stream objects we've been using
      delete mMidFluidCollision;
      delete mWallCollision;
      delete mInletCollision;
      delete mOutletCollision;
      delete mInletWallCollision;
      delete mOutletWallCollision;

      // Delete various other arrays used
      delete[] inlet_count;
      delete[] inlet_normal;
      delete[] average_inlet_velocity;
      delete[] peak_inlet_velocity;
    }
  }
}
