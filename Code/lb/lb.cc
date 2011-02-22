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

    // Set up of min/max values at the beginning of each pulsatile cycle.
    void LBM::InitMinMaxValues(void)
    {
      // All values will be in [0, infinity) so it's adequate
      // to set initial maxima to -1, and minima to std::numeric_limits<double>::max();
      mMinsAndMaxes.MaxDensity = -1.0;
      mMinsAndMaxes.MaxVelocity = -1.0;
      mMinsAndMaxes.MaxStress = -1.0;

      mMinsAndMaxes.MinDensity = std::numeric_limits<double>::max();
      mMinsAndMaxes.MinVelocity = std::numeric_limits<double>::max();
      mMinsAndMaxes.MinStress = std::numeric_limits<double>::max();
    }

    // Calculate the BCs for each boundary site type and the
    // non-equilibrium distribution functions.
    void LBM::CalculateBC(double f[],
                          hemelb::lb::SiteType iSiteType,
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

      if (iSiteType == hemelb::lb::FLUID_TYPE)
      {
        D3Q15::CalculateDensityAndVelocity(f, *density, *vx, *vy, *vz);
      }
      else
      {
        if (iSiteType == hemelb::lb::INLET_TYPE)
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

    const hemelb::lb::LbmParameters *LBM::GetLbmParams()
    {
      return &mParams;
    }

    LBM::LBM(hemelb::SimConfig *iSimulationConfig,
             const hemelb::topology::NetworkTopology * iNetTop,
             hemelb::lb::GlobalLatticeData &bGlobLatDat,
             int iSteeringSessionId,
             double * oFileReadTime)
    {
      steering_session_id = iSteeringSessionId;
      period = iSimulationConfig->StepsPerCycle;
      voxel_size = iSimulationConfig->VoxelSize;

      mNetTopology = iNetTop;
      mSimConfig = iSimulationConfig;

      double fileReadStartTime = hemelb::util::myClock();
      ReadConfig(bGlobLatDat);
      *oFileReadTime = hemelb::util::myClock() - fileReadStartTime;

      ReadParameters();

      InitCollisions();
    }

    void LBM::CalculateMouseFlowField(hemelb::vis::ColPixel *col_pixel_p,
                                      double &mouse_pressure,
                                      double &mouse_stress,
                                      double density_threshold_min,
                                      double density_threshold_minmax_inv,
                                      double stress_threshold_max_inv)
    {
      double density = density_threshold_min + col_pixel_p->density / density_threshold_minmax_inv;
      double stress = col_pixel_p->stress / stress_threshold_max_inv;

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

    void LBM::Initialise(int* iFTranslator)
    {
      receivedFTranslator = iFTranslator;
    }

    void LBM::SetInitialConditions(hemelb::lb::LocalLatticeData &bLocalLatDat)
    {
      double *f_old_p, *f_new_p, f_eq[D3Q15::NUMVECTORS];
      double density;

      density = 0.;

      for (int i = 0; i < outlets; i++)
      {
        density += outlet_density_avg[i] - outlet_density_amp[i];
      }
      density /= outlets;

      for (int i = 0; i < bLocalLatDat.GetLocalFluidSiteCount(); i++)
      {
        D3Q15::CalculateFeq(density, 0.0, 0.0, 0.0, f_eq);

        f_old_p = &bLocalLatDat.FOld[i * D3Q15::NUMVECTORS];
        f_new_p = &bLocalLatDat.FNew[i * D3Q15::NUMVECTORS];

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

    void LBM::RequestComms(net::Net* net, lb::LocalLatticeData* bLocalLatDat)
    {
      for (std::vector<hemelb::topology::NeighbouringProcessor*>::const_iterator it =
          mNetTopology->NeighbouringProcs.begin(); it != mNetTopology->NeighbouringProcs.end(); it++)
      {
        // Request the receive into the appropriate bit of FOld.
        net->RequestReceive<double> (&bLocalLatDat->FOld[ (*it)->FirstSharedF],
                                      (*it)->SharedFCount, (*it)->Rank);

        // Request the send from the right bit of
        net->RequestSend<double> (&bLocalLatDat->FNew[ (*it)->FirstSharedF], (*it)->SharedFCount,
                                   (*it)->Rank);
      }
    }

    void LBM::PreSend(lb::LocalLatticeData* bLocalLatDat, int perform_rt)
    {
      int offset = bLocalLatDat->my_inner_sites;

      for (int collision_type = 0; collision_type < COLLISION_TYPES; collision_type++)
      {
        GetCollision(collision_type)->DoCollisions(
                                                   perform_rt,
                                                   offset,
                                                   bLocalLatDat->my_inter_collisions[collision_type],
                                                   mParams, mMinsAndMaxes, *bLocalLatDat,
                                                   hemelb::vis::controller);
        offset += bLocalLatDat->my_inter_collisions[collision_type];
      }
    }

    void LBM::PreReceive(int perform_rt, lb::LocalLatticeData* bLocalLatDat)
    {
      int offset = 0;

      for (int collision_type = 0; collision_type < COLLISION_TYPES; collision_type++)
      {
        GetCollision(collision_type)->DoCollisions(
                                                   perform_rt,
                                                   offset,
                                                   bLocalLatDat->my_inner_collisions[collision_type],
                                                   mParams, mMinsAndMaxes, *bLocalLatDat,
                                                   hemelb::vis::controller);
        offset += bLocalLatDat->my_inner_collisions[collision_type];
      }
    }

    void LBM::PostReceive(lb::LocalLatticeData* bLocalLatDat, int perform_rt)
    {
      // Copy the distribution functions received from the neighbouring
      // processors into the destination buffer "f_new".
      for (int i = 0; i < mNetTopology->TotalSharedFs; i++)
      {
        bLocalLatDat->FNew[receivedFTranslator[i]]
            = bLocalLatDat->FOld[mNetTopology->NeighbouringProcs[0]->FirstSharedF + i];
      }

      // Do any cleanup steps necessary on boundary nodes
      int offset = 0;

      for (int collision_type = 0; collision_type < COLLISION_TYPES; collision_type++)
      {
        GetCollision(collision_type)->PostStep(perform_rt, offset,
                                               bLocalLatDat->my_inner_collisions[collision_type],
                                               mParams, mMinsAndMaxes, *bLocalLatDat,
                                               hemelb::vis::controller);
        offset += bLocalLatDat->my_inner_collisions[collision_type];
      }

      for (int collision_type = 0; collision_type < COLLISION_TYPES; collision_type++)
      {
        GetCollision(collision_type)->PostStep(perform_rt, offset,
                                               bLocalLatDat->my_inter_collisions[collision_type],
                                               mParams, mMinsAndMaxes, *bLocalLatDat,
                                               hemelb::vis::controller);
        offset += bLocalLatDat->my_inter_collisions[collision_type];
      }
    }

    void LBM::EndIteration(lb::LocalLatticeData* bLocalLatDat)
    {
      // Swap f_old and f_new ready for the next timestep.
      double *temp = bLocalLatDat->FOld;
      bLocalLatDat->FOld = bLocalLatDat->FNew;
      bLocalLatDat->FNew = temp;
    }

    void LBM::CalculateFlowFieldValues()
    {
      double *local_data;
      double *global_data;

      int i;

      int lMaxInlets = hemelb::util::max(6 + inlets, 2 * inlets);

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
      MPI_Reduce(&local_data[3], &global_data[3], 3 + inlets, MPI_DOUBLE, MPI_MAX, 0,
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

      period = period;

      inlets = inlets;
    }

    // Update peak and average inlet velocities local to the current subdomain.
    void LBM::UpdateInletVelocities(int time_step,
                                    lb::LocalLatticeData &iLocalLatDat,
                                    net::Net *net)
    {
      double density;
      double vx, vy, vz;
      double velocity;

      int offset;
      int inlet_id;
      int i;
      int c1, c2;

      if (time_step == 1)
      {
        for (i = 0; i < inlets; i++)
        {
          peak_inlet_velocity[i] = -BIG_NUMBER;
          average_inlet_velocity[i] = 0.;
          inlet_count[i] = 0;
        }
      }

      c1 = 15;
      c2 = 0;

      offset = iLocalLatDat.my_inner_collisions[0] + iLocalLatDat.my_inner_collisions[1];

      for (i = offset; i < offset + iLocalLatDat.my_inner_collisions[2]; i++)
      {
        D3Q15::CalculateDensityAndVelocity(&iLocalLatDat.FOld[i * c1 + c2], density, vx, vy, vz);

        inlet_id = iLocalLatDat.GetBoundaryId(i);

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
      offset = iLocalLatDat.my_inner_sites + iLocalLatDat.my_inter_collisions[0]
          + iLocalLatDat.my_inter_collisions[1];

      for (i = offset; i < offset + iLocalLatDat.my_inter_collisions[2]; i++)
      {
        D3Q15::CalculateDensityAndVelocity(&iLocalLatDat.FOld[i * c1 + c2], density, vx, vy, vz);

        inlet_id = iLocalLatDat.GetBoundaryId(i);

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
    void LBM::Restart(hemelb::lb::LocalLatticeData &iLocalLatDat)
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

      SetInitialConditions(iLocalLatDat);
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
