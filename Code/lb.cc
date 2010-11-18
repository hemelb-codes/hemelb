// In this file, the functions useful to calculate the equilibrium distribution
// function, momentums, the effective von Mises stress and the boundary conditions
// are reported

#include <math.h>
#include <limits>

#include "lb.h"
#include "utilityFunctions.h"
#include "vis/RayTracer.h"

double LBM::lbmConvertPressureToLatticeUnits(double pressure)
{
  return Cs2 + (pressure - REFERENCE_PRESSURE) * mmHg_TO_PASCAL
      * (PULSATILE_PERIOD / (period * voxel_size)) * (PULSATILE_PERIOD
      / (period * voxel_size)) / BLOOD_DENSITY;
}

double LBM::lbmConvertPressureToPhysicalUnits(double pressure)
{
  return REFERENCE_PRESSURE + ( (pressure / Cs2 - 1.0) * Cs2) * BLOOD_DENSITY
      * ( (period * voxel_size) / PULSATILE_PERIOD) * ( (period * voxel_size)
      / PULSATILE_PERIOD) / mmHg_TO_PASCAL;
}

double LBM::lbmConvertPressureGradToLatticeUnits(double pressure_grad)
{
  return pressure_grad * mmHg_TO_PASCAL * (PULSATILE_PERIOD / (period
      * voxel_size)) * (PULSATILE_PERIOD / (period * voxel_size))
      / BLOOD_DENSITY;
}

double LBM::lbmConvertPressureGradToPhysicalUnits(double pressure_grad)
{
  return pressure_grad * BLOOD_DENSITY * ( (period * voxel_size)
      / PULSATILE_PERIOD) * ( (period * voxel_size) / PULSATILE_PERIOD)
      / mmHg_TO_PASCAL;
}

double LBM::lbmConvertVelocityToLatticeUnits(double velocity)
{
  return velocity * ( ( (mParams.Tau - 0.5) / 3.0) * voxel_size)
      / (BLOOD_VISCOSITY / BLOOD_DENSITY);
}

double LBM::lbmConvertVelocityToPhysicalUnits(double velocity)
{
  // convert velocity from lattice units to physical units (m/s)
  return velocity * (BLOOD_VISCOSITY / BLOOD_DENSITY) / ( ( (mParams.Tau - 0.5)
      / 3.0) * voxel_size);
}

double LBM::lbmConvertStressToLatticeUnits(double stress)
{
  return stress * (BLOOD_DENSITY / (BLOOD_VISCOSITY * BLOOD_VISCOSITY))
      * ( ( (mParams.Tau - 0.5) / 3.0) * voxel_size) * ( ( (mParams.Tau - 0.5)
      / 3.0) * voxel_size);
}

double LBM::lbmConvertStressToPhysicalUnits(double stress)
{
  // convert stress from lattice units to physical units (Pa)
  return stress * BLOOD_VISCOSITY * BLOOD_VISCOSITY / (BLOOD_DENSITY
      * ( ( (mParams.Tau - 0.5) / 3.0) * voxel_size) * ( ( (mParams.Tau - 0.5)
      / 3.0) * voxel_size));
}

void LBM::RecalculateTauViscosityOmega()
{
  mParams.Tau = 0.5 + (PULSATILE_PERIOD * BLOOD_VISCOSITY / BLOOD_DENSITY)
      / (Cs2 * period * voxel_size * voxel_size);

  mParams.Omega = -1.0 / mParams.Tau;
  mParams.StressParameter = (1.0 - 1.0 / (2.0 * mParams.Tau)) / sqrt(2.0);
}

// Set up of min/max values at the beginning of each pulsatile cycle.
void LBM::lbmInitMinMaxValues(void)
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
void LBM::lbmCalculateBC(double f[],
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

void LBM::lbmUpdateBoundaryDensities(int cycle_id, int time_step)
{
  double w = 2.0 * PI / period;

  for (int i = 0; i < inlets; i++)
  {
    inlet_density[i] = inlet_density_avg[i] + inlet_density_amp[i] * cos(w
        * (double) time_step + inlet_density_phs[i]);
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

void LBM::lbmInit(hemelb::SimConfig *iSimulationConfig,
                  const hemelb::topology::NetworkTopology * iNetTop,
                  hemelb::lb::GlobalLatticeData &bGlobLatDat,
                  int iSteeringSessionId,
                  int iPeriod,
                  double iVoxelSize,
                  Net *net)
{
  steering_session_id = iSteeringSessionId;
  period = iPeriod;
  voxel_size = iVoxelSize;

  mNetTopology = iNetTop;
  mSimConfig = iSimulationConfig;

  lbmReadConfig(net, bGlobLatDat);

  lbmReadParameters();

  lbmInitCollisions();
}

void LBM::CalculateMouseFlowField(hemelb::vis::ColPixel *col_pixel_p,
                                  double &mouse_pressure,
                                  double &mouse_stress,
                                  double density_threshold_min,
                                  double density_threshold_minmax_inv,
                                  double stress_threshold_max_inv)
{
  double density = density_threshold_min + col_pixel_p->density
      / density_threshold_minmax_inv;
  double stress = col_pixel_p->stress / stress_threshold_max_inv;

  mouse_pressure = lbmConvertPressureToPhysicalUnits(density * Cs2);
  mouse_stress = lbmConvertStressToPhysicalUnits(stress);
}

void LBM::lbmInitCollisions()
{
  // TODO Note that the convergence checking is not yet implemented in the
  // new boundary condition hierarchy system.
  // It'd be nice to do this with something like
  // MidFluidCollision = new ConvergenceCheckingWrapper(new WhateverMidFluidCollision());

  mMidFluidCollision = new hemelb::lb::collisions::ImplSimpleCollideAndStream();
  mWallCollision = new hemelb::lb::collisions::ImplZeroVelocityEquilibrium();
  mInletCollision
      = new hemelb::lb::collisions::ImplNonZeroVelocityBoundaryDensity(
                                                                       inlet_density);
  mOutletCollision
      = new hemelb::lb::collisions::ImplNonZeroVelocityBoundaryDensity(
                                                                       outlet_density);
  mInletWallCollision
      = new hemelb::lb::collisions::ImplZeroVelocityBoundaryDensity(
                                                                    inlet_density);
  mOutletWallCollision
      = new hemelb::lb::collisions::ImplZeroVelocityBoundaryDensity(
                                                                    outlet_density);
}

void LBM::lbmSetInitialConditions(hemelb::lb::LocalLatticeData &bLocalLatDat)
{
  double *f_old_p, *f_new_p, f_eq[D3Q15::NUMVECTORS];
  double density;

  density = 0.;

  for (int i = 0; i < outlets; i++)
  {
    density += outlet_density_avg[i] - outlet_density_amp[i];
  }
  density /= outlets;

  for (int i = 0; i < bLocalLatDat.LocalFluidSites; i++)
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

// The entire simulation time step takes place through this function
// when the convergence criterion is not applied. Communications
// automatically handle the streaming stage pertaining to neighbouring
// subdomains.
hemelb::lb::Stability LBM::lbmCycle(int perform_rt,
                  Net *net,
                  hemelb::lb::LocalLatticeData &bLocalLatDat,
                  double &bLbTime,
                  double &bMPISendTime,
                  double &bMPIWaitTime)
{
  net->ReceiveFromNeighbouringProcessors(bLocalLatDat);

  int offset = net->my_inner_sites;

#ifndef NOMPI
  double lPreLbTimeOne = MPI_Wtime();
#endif

  for (int collision_type = 0; collision_type < COLLISION_TYPES; collision_type++)
  {
    GetCollision(collision_type)->DoCollisions(
                                               perform_rt,
                                               offset,
                                               net->my_inter_collisions[collision_type],
                                               mParams, mMinsAndMaxes,
                                               bLocalLatDat,
                                               hemelb::vis::controller);
    offset += net->my_inter_collisions[collision_type];
  }

#ifndef NOMPI
  double lPreSendTime = MPI_Wtime();
  bLbTime += (lPreSendTime - lPreLbTimeOne);
#endif

  net->SendToNeighbouringProcessors(bLocalLatDat);

#ifndef NOMPI
  double lPreLbTimeTwo = MPI_Wtime();
  bMPISendTime += (lPreLbTimeTwo - lPreSendTime);
#endif

  offset = 0;

  for (int collision_type = 0; collision_type < COLLISION_TYPES; collision_type++)
  {
    GetCollision(collision_type)->DoCollisions(
                                               perform_rt,
                                               offset,
                                               net->my_inner_collisions[collision_type],
                                               mParams, mMinsAndMaxes,
                                               bLocalLatDat,
                                               hemelb::vis::controller);
    offset += net->my_inner_collisions[collision_type];
  }

#ifndef NOMPI
  double lPreMPIWaitTime = MPI_Wtime();
  bLbTime += (lPreMPIWaitTime - lPreLbTimeTwo);
#endif

  net->UseDataFromNeighbouringProcs(bLocalLatDat);

#ifndef NOMPI
  double lPrePostStepTime = MPI_Wtime();
  bMPIWaitTime += (lPrePostStepTime - lPreMPIWaitTime);
#endif

  // Do any cleanup steps necessary on boundary nodes
  offset = 0;

  for (int collision_type = 0; collision_type < COLLISION_TYPES; collision_type++)
  {
    GetCollision(collision_type)->PostStep(
                                           perform_rt,
                                           offset,
                                           net->my_inner_collisions[collision_type],
                                           mParams, mMinsAndMaxes,
                                           bLocalLatDat,
                                           hemelb::vis::controller);
    offset += net->my_inner_collisions[collision_type];
  }

  for (int collision_type = 0; collision_type < COLLISION_TYPES; collision_type++)
  {
    GetCollision(collision_type)->PostStep(
                                           perform_rt,
                                           offset,
                                           net->my_inter_collisions[collision_type],
                                           mParams, mMinsAndMaxes,
                                           bLocalLatDat,
                                           hemelb::vis::controller);
    offset += net->my_inter_collisions[collision_type];
  }

#ifndef NOMPI
  bLbTime += (MPI_Wtime() - lPrePostStepTime);
#endif

  // Swap f_old and f_new ready for the next timestep.
  double *temp = bLocalLatDat.FOld;
  bLocalLatDat.FOld = bLocalLatDat.FNew;
  bLocalLatDat.FNew = temp;

  return hemelb::lb::Stable;
}

void LBM::lbmCalculateFlowFieldValues()
{
  double *local_data;
  double *global_data;

  int i;

  int lMaxInlets = hemelb::util::max(6 + inlets, 2 * inlets);

  local_data = new double[lMaxInlets];
  global_data = new double[lMaxInlets];

#ifndef NOMPI
  local_data[0] = mMinsAndMaxes.MinDensity;
  local_data[1] = mMinsAndMaxes.MinVelocity;
  local_data[2] = mMinsAndMaxes.MinStress;
  local_data[3] = mMinsAndMaxes.MaxDensity;
  local_data[4] = mMinsAndMaxes.MaxVelocity;
  local_data[5] = mMinsAndMaxes.MaxStress;

  memcpy(&local_data[6], lbm_peak_inlet_velocity, sizeof(double) * inlets);

  MPI_Reduce(&local_data[0], &global_data[0], 3, MPI_DOUBLE, MPI_MIN, 0,
             MPI_COMM_WORLD);
  MPI_Reduce(&local_data[3], &global_data[3], 3 + inlets, MPI_DOUBLE, MPI_MAX,
             0, MPI_COMM_WORLD);

  mMinsAndMaxes.MinDensity = global_data[0];
  mMinsAndMaxes.MinVelocity = global_data[1];
  mMinsAndMaxes.MinStress = global_data[2];
  mMinsAndMaxes.MaxDensity = global_data[3];
  mMinsAndMaxes.MaxVelocity = global_data[4];
  mMinsAndMaxes.MaxStress = global_data[5];

  memcpy(lbm_peak_inlet_velocity, &global_data[6], sizeof(double) * inlets);

  for (i = 0; i < inlets; i++)
  {
    local_data[i] = lbm_average_inlet_velocity[i];
    local_data[inlets + i] = lbm_inlet_count[i];
  }
  MPI_Reduce(local_data, global_data, 2 * inlets, MPI_DOUBLE, MPI_SUM, 0,
             MPI_COMM_WORLD);

  for (i = 0; i < inlets; i++)
  {
    lbm_average_inlet_velocity[i] = global_data[i];
    lbm_inlet_count[i] = global_data[inlets + i];
  }
#endif // NOMPI
  delete[] global_data;
  delete[] local_data;

  for (i = 0; i < inlets; i++)
  {
    lbm_average_inlet_velocity[i] /= lbm_inlet_count[i];
    lbm_average_inlet_velocity[i]
        = lbmConvertVelocityToPhysicalUnits(lbm_average_inlet_velocity[i]);
    lbm_peak_inlet_velocity[i]
        = lbmConvertVelocityToPhysicalUnits(lbm_peak_inlet_velocity[i]);
  }

  period = period;

  inlets = inlets;
}

int LBM::IsUnstable(hemelb::lb::LocalLatticeData &iLocalLatDat, Net *net)
{
  int is_unstable, stability;

  is_unstable = 0;

  for (int i = 0; i < iLocalLatDat.LocalFluidSites; i++)
  {
    for (unsigned int l = 0; l < D3Q15::NUMVECTORS; l++)
    {
      if (iLocalLatDat.FOld[i * D3Q15::NUMVECTORS + l] < 0.)
      {
        is_unstable = 1;
      }
    }
  }

#ifndef NOMPI
  net->err = MPI_Allreduce(&is_unstable, &stability, 1, MPI_INT, MPI_MAX,
                           MPI_COMM_WORLD);
  is_unstable = stability;
#endif

  return is_unstable;
}

// Update peak and average inlet velocities local to the current subdomain. 
void LBM::lbmUpdateInletVelocities(int time_step,
                                   hemelb::lb::LocalLatticeData &iLocalLatDat,
                                   Net *net)
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
      lbm_peak_inlet_velocity[i] = -BIG_NUMBER;
      lbm_average_inlet_velocity[i] = 0.;
      lbm_inlet_count[i] = 0;
    }
  }

  c1 = 15;
  c2 = 0;

  offset = net->my_inner_collisions[0] + net->my_inner_collisions[1];

  for (i = offset; i < offset + net->my_inner_collisions[2]; i++)
  {
    D3Q15::CalculateDensityAndVelocity(&iLocalLatDat.FOld[i * c1 + c2],
                                       density, vx, vy, vz);

    inlet_id = iLocalLatDat.GetBoundaryId(i);

    if (is_inlet_normal_available)
    {
      vx *= lbm_inlet_normal[3 * inlet_id + 0];
      vy *= lbm_inlet_normal[3 * inlet_id + 1];
      vz *= lbm_inlet_normal[3 * inlet_id + 2];

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
    lbm_peak_inlet_velocity[inlet_id] = fmax(lbm_peak_inlet_velocity[inlet_id],
                                             velocity);
    lbm_average_inlet_velocity[inlet_id] += velocity;
    ++lbm_inlet_count[inlet_id];
  }
  offset = net->my_inner_sites + net->my_inter_collisions[0]
      + net->my_inter_collisions[1];

  for (i = offset; i < offset + net->my_inter_collisions[2]; i++)
  {
    D3Q15::CalculateDensityAndVelocity(&iLocalLatDat.FOld[i * c1 + c2],
                                       density, vx, vy, vz);

    inlet_id = iLocalLatDat.GetBoundaryId(i);

    if (is_inlet_normal_available)
    {
      vx *= lbm_inlet_normal[3 * inlet_id + 0];
      vy *= lbm_inlet_normal[3 * inlet_id + 1];
      vz *= lbm_inlet_normal[3 * inlet_id + 2];

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
    lbm_peak_inlet_velocity[inlet_id] = fmax(lbm_peak_inlet_velocity[inlet_id],
                                             velocity);
    lbm_average_inlet_velocity[inlet_id] += velocity;
    ++lbm_inlet_count[inlet_id];
  }
}

double LBM::GetAverageInletVelocity(int iInletNumber)
{
  return lbm_average_inlet_velocity[iInletNumber];
}
double LBM::GetPeakInletVelocity(int iInletNumber)
{
  return lbm_peak_inlet_velocity[iInletNumber];
}

// In the case of instability, this function restart the simulation
// with twice as many time steps per period and update the parameters
// that depends on this change.
void LBM::lbmRestart(hemelb::lb::LocalLatticeData &iLocalLatDat)
{
  int i;

  for (i = 0; i < inlets; i++)
  {
    inlet_density_avg[i]
        = lbmConvertPressureToPhysicalUnits(inlet_density_avg[i] * Cs2);
    inlet_density_amp[i]
        = lbmConvertPressureGradToPhysicalUnits(inlet_density_amp[i] * Cs2);
  }
  for (i = 0; i < outlets; i++)
  {
    outlet_density_avg[i]
        = lbmConvertPressureToPhysicalUnits(outlet_density_avg[i] * Cs2);
    outlet_density_amp[i]
        = lbmConvertPressureGradToPhysicalUnits(outlet_density_amp[i] * Cs2);
  }
  period *= 2;

  for (i = 0; i < inlets; i++)
  {
    inlet_density_avg[i]
        = lbmConvertPressureToLatticeUnits(inlet_density_avg[i]) / Cs2;
    inlet_density_amp[i]
        = lbmConvertPressureGradToLatticeUnits(inlet_density_amp[i]) / Cs2;
  }
  for (i = 0; i < outlets; i++)
  {
    outlet_density_avg[i]
        = lbmConvertPressureToLatticeUnits(outlet_density_avg[i]) / Cs2;
    outlet_density_amp[i]
        = lbmConvertPressureGradToLatticeUnits(outlet_density_amp[i]) / Cs2;
  }

  RecalculateTauViscosityOmega();

  lbmSetInitialConditions(iLocalLatDat);
}

double LBM::GetMinPhysicalPressure()
{
  return lbmConvertPressureToPhysicalUnits(mMinsAndMaxes.MinDensity * Cs2);
}
double LBM::GetMaxPhysicalPressure()
{
  return lbmConvertPressureToPhysicalUnits(mMinsAndMaxes.MaxDensity * Cs2);
}
double LBM::GetMinPhysicalVelocity()
{
  return lbmConvertVelocityToPhysicalUnits(mMinsAndMaxes.MinVelocity);
}
double LBM::GetMaxPhysicalVelocity()
{
  return lbmConvertVelocityToPhysicalUnits(mMinsAndMaxes.MaxVelocity);
}
double LBM::GetMinPhysicalStress()
{
  return lbmConvertStressToPhysicalUnits(mMinsAndMaxes.MinStress);
}
double LBM::GetMaxPhysicalStress()
{
  return lbmConvertStressToPhysicalUnits(mMinsAndMaxes.MaxStress);
}

LBM::~LBM()
{
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
  delete[] lbm_inlet_count;
  delete[] lbm_inlet_normal;
  delete[] lbm_average_inlet_velocity;
  delete[] lbm_peak_inlet_velocity;
}
