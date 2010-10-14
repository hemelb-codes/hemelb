// In this file, the functions useful to calculate the equilibrium distribution
// function, momentums, the effective von Mises stress and the boundary conditions
// are reported

#include <math.h>

#include "lb.h"
#include "utilityFunctions.h"
#include "vis/RayTracer.h"

void (*lbmUpdateSiteData[2])(bool isPostStep, double omega, int i, double *density, double *vx,
  double *vy, double *vz, double *velocity, Net* net, LBM* iLbm,
  hemelb::lb::collisions::Collision* iCollision);

double LBM::lbmConvertPressureToLatticeUnits (double pressure)
{
  return Cs2 + (pressure - REFERENCE_PRESSURE) * mmHg_TO_PASCAL *
    (PULSATILE_PERIOD / (period * voxel_size)) *
    (PULSATILE_PERIOD / (period * voxel_size)) / BLOOD_DENSITY;
}


double LBM::lbmConvertPressureToPhysicalUnits (double pressure)
{
  return REFERENCE_PRESSURE + ((pressure / Cs2 - 1.0) * Cs2) * BLOOD_DENSITY *
    ((period * voxel_size) / PULSATILE_PERIOD) *
    ((period * voxel_size) / PULSATILE_PERIOD) / mmHg_TO_PASCAL;
}


double LBM::lbmConvertPressureGradToLatticeUnits (double pressure_grad)
{
  return pressure_grad * mmHg_TO_PASCAL *
    (PULSATILE_PERIOD / (period * voxel_size)) *
    (PULSATILE_PERIOD / (period * voxel_size)) / BLOOD_DENSITY;
}


double LBM::lbmConvertPressureGradToPhysicalUnits (double pressure_grad)
{
  return pressure_grad * BLOOD_DENSITY *
    ((period * voxel_size) / PULSATILE_PERIOD) *
    ((period * voxel_size) / PULSATILE_PERIOD) / mmHg_TO_PASCAL;
}


double LBM::lbmConvertVelocityToLatticeUnits (double velocity)
{
  return velocity * (((tau-0.5)/3.0) * voxel_size) / (BLOOD_VISCOSITY / BLOOD_DENSITY);
}


double LBM::lbmConvertVelocityToPhysicalUnits (double velocity)
{
  // convert velocity from lattice units to physical units (m/s)
  return velocity * (BLOOD_VISCOSITY / BLOOD_DENSITY) / (((tau-0.5)/3.0) * voxel_size);
}


double LBM::lbmConvertStressToLatticeUnits (double stress)
{
  return stress * (BLOOD_DENSITY / (BLOOD_VISCOSITY * BLOOD_VISCOSITY)) *
    (((tau-0.5)/3.0) * voxel_size) *
    (((tau-0.5)/3.0) * voxel_size);
}


double LBM::lbmConvertStressToPhysicalUnits (double stress)
{
  // convert stress from lattice units to physical units (Pa)
  return stress * BLOOD_VISCOSITY * BLOOD_VISCOSITY /
    (BLOOD_DENSITY * (((tau-0.5)/3.0) * voxel_size) * (((tau-0.5)/3.0) * voxel_size));
}


void LBM::RecalculateTauViscosityOmega ()
{
  tau = 0.5 + (PULSATILE_PERIOD * BLOOD_VISCOSITY / BLOOD_DENSITY) /
    (Cs2 * period * voxel_size * voxel_size);

  viscosity = ((2.0 * tau - 1.0) / 6.0);
  omega = -1.0 / tau;
  lbm_stress_par = (1.0 - 1.0 / (2.0 * tau)) / sqrt(2.0);
}

// Set up of min/max values at the beginning of each pulsatile cycle.
void LBM::lbmInitMinMaxValues (void)
{
  lbm_density_min = +1.0e+30;
  lbm_density_max = +1.0e-30;

  lbm_velocity_min = +1.0e+30;
  lbm_velocity_max = +1.0e-30;

  lbm_stress_min = +1.0e+30;
  lbm_stress_max = +1.0e-30;
}

// Update the min/max values local to the current subdomain.
void LBM::lbmUpdateMinMaxValues (double density, double velocity, double stress)
{
  lbm_density_min = (density < lbm_density_min) ? density : lbm_density_min;
  lbm_density_max = (density > lbm_density_max) ? density : lbm_density_max;

  lbm_velocity_min = (velocity < lbm_velocity_min) ? velocity : lbm_velocity_min;
  lbm_velocity_max = (velocity > lbm_velocity_max) ? velocity : lbm_velocity_max;

  lbm_stress_min = (stress < lbm_stress_min) ? stress : lbm_stress_min;
  lbm_stress_max = (stress > lbm_stress_max) ? stress : lbm_stress_max;
}

// Fluid site updating for full-production runs.
void lbmUpdateSiteDataSim (const bool isPostStep, double omega, int i, double *density, double *vx,double *vy, double *vz, double *velocity, Net* net,
  LBM* iLbm, hemelb::lb::collisions::Collision* iCollision)
{
  double f_neq[D3Q15::NUMVECTORS];
  double stress;

  bool lUpdateMinsAndMaxes = true;

  if(isPostStep)
  {
    lUpdateMinsAndMaxes = iCollision->PostStep(omega, i, density, vx, vy, vz, f_neq, net);
  }
  else
  {
    iCollision->DoCollisions(omega, i, density, vx, vy, vz, f_neq, net);
  }

  if (lUpdateMinsAndMaxes)
  {

    *vx *= (1.0 / *density);
    *vy *= (1.0 / *density);
    *vz *= (1.0 / *density);
    *velocity = sqrt(*vx * *vx + *vy * *vy + *vz * *vz);

    if (lbm_stress_type == SHEAR_STRESS)
    {
      if (net->GetNormalToWall(i)[0] > 1.0e+30)
      {
        stress = 0.0;
      }
      else
      {
        D3Q15::CalculateShearStress(*density, f_neq, net->GetNormalToWall(i), &stress,
                                    lbm_stress_par);
      }
    }
    else
    {
      D3Q15::CalculateVonMisesStress(f_neq, &stress, lbm_stress_par);
    }
    iLbm->lbmUpdateMinMaxValues(*density, *velocity, stress);
  }
}

// Fluid site updating for full-production runs plus computation of flow field values for visualisation purposes.
void lbmUpdateSiteDataSimPlusVis (const bool isPostStep, double omega, int i, double *density, double *vx,double *vy, double *vz, double *velocity, Net* net,
  LBM *iLbm, hemelb::lb::collisions::Collision* iCollision)
{
  double f_neq[D3Q15::NUMVECTORS];
  double stress;

  bool lUpdateMinsAndMaxes = true;

  if(isPostStep)
  {
    lUpdateMinsAndMaxes = iCollision->PostStep(omega, i, density, vx, vy, vz, f_neq, net);
  }
  else
  {
    iCollision->DoCollisions(omega, i, density, vx, vy, vz, f_neq, net);
  }

  if (lUpdateMinsAndMaxes)
  {
    *vx *= (1.0 / *density);
    *vy *= (1.0 / *density);
    *vz *= (1.0 / *density);
    *velocity = sqrt(*vx * *vx + *vy * *vy + *vz * *vz);

    if (lbm_stress_type == SHEAR_STRESS)
    {
      if (net->GetNormalToWall(i)[0] >= 1.0e+30)
      {
        iLbm->lbmUpdateMinMaxValues(*density, *velocity, 0.0);
        hemelb::vis::rtUpdateClusterVoxel(i, *density, *velocity, 1.0e+30F);
      }
      else
      {
        D3Q15::CalculateShearStress(*density, f_neq, net->GetNormalToWall(i), &stress,
                                    lbm_stress_par);
        iLbm->lbmUpdateMinMaxValues(*density, *velocity, stress);
        hemelb::vis::rtUpdateClusterVoxel(i, *density, *velocity, stress);
      }
    }
    else
    {
      D3Q15::CalculateVonMisesStress(f_neq, &stress, lbm_stress_par);
      iLbm->lbmUpdateMinMaxValues(*density, *velocity, stress);
      hemelb::vis::rtUpdateClusterVoxel(i, *density, *velocity, stress);
    }
  }
}

// Calculate the BCs for each boundary site type and the
// non-equilibrium distribution functions.
void LBM::lbmCalculateBC (double f[], unsigned int site_data, double *density,
		     double *vx, double *vy, double *vz, double f_neq[])
{
  double dummy_density;

  int l;

  unsigned int boundary_type, boundary_id;

  for (l = 0; l < D3Q15::NUMVECTORS; l++)
    {
      f_neq[ l ] = f[ l ];
    }
  boundary_type = site_data & SITE_TYPE_MASK;

  if (boundary_type == FLUID_TYPE)
    {
      D3Q15::CalculateDensityAndVelocity(f, density, vx, vy, vz);
    }
  else
    {
      boundary_id     = (site_data & BOUNDARY_ID_MASK)     >> BOUNDARY_ID_SHIFT;

      if (boundary_type == INLET_TYPE)
	{
	  *density = inlet_density[ boundary_id ];
	}
      else
	{
	  *density = outlet_density[ boundary_id ];
	}

      D3Q15::CalculateDensityAndVelocity(f, &dummy_density, vx, vy, vz);
      D3Q15::CalculateFeq(*density, *vx, *vy, *vz, f);

    }
  for (l = 0; l < D3Q15::NUMVECTORS; l++)
    {
      f_neq[ l ] -= f[ l ];
    }
}


void LBM::lbmUpdateBoundaryDensities (int cycle_id, int time_step)
{
  double w = 2.0 * PI / period;

  for (int i = 0; i < inlets; i++)
    {
      inlet_density[i] = inlet_density_avg[i] + inlet_density_amp[i] * cos(w * (double)time_step + inlet_density_phs[i]);
    }
  for (int i = 0; i < outlets; i++)
    {
      outlet_density[i] = outlet_density_avg[i] + outlet_density_amp[i] * cos(w * (double)time_step + outlet_density_phs[i]);
    }
}


void LBM::lbmInit (char *system_file_name_in, char *parameters_file_name, Net *net)
{
  system_file_name = system_file_name_in;

  lbmUpdateSiteData[0] = lbmUpdateSiteDataSim;
  lbmUpdateSiteData[1] = lbmUpdateSiteDataSimPlusVis;

  lbm_terminate_simulation = 0;

  lbmReadConfig (net);

  lbmReadParameters(parameters_file_name, net);

  lbmInitCollisions();
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
      = new hemelb::lb::collisions::ImplNonZeroVelocityBoundaryDensity(inlet_density);
  mOutletCollision
      = new hemelb::lb::collisions::ImplNonZeroVelocityBoundaryDensity(outlet_density);
  mInletWallCollision
      = new hemelb::lb::collisions::ImplZeroVelocityBoundaryDensity(inlet_density);
  mOutletWallCollision
      = new hemelb::lb::collisions::ImplZeroVelocityBoundaryDensity(outlet_density);
}

void LBM::lbmSetInitialConditions (Net *net)
{
  double *f_old_p, *f_new_p, f_eq[D3Q15::NUMVECTORS];
  double density;
  int i, l;

  density = 0.;

  for (i = 0; i < outlets; i++)
    {
      density += outlet_density_avg[i] - outlet_density_amp[i];
    }
  density /= outlets;

  for (i = 0; i < net->my_sites; i++)
    {
    D3Q15::CalculateFeq(density, 0.0, 0.0, 0.0, f_eq);

	  f_old_p = &f_old[ i*D3Q15::NUMVECTORS ];
	  f_new_p = &f_new[ i*D3Q15::NUMVECTORS ];

	  for (l = 0; l < D3Q15::NUMVECTORS; l++)
	    {
	      f_new_p[ l ] = f_old_p[ l ] = f_eq[ l ];
	    }
    }
}

// TODO HACK
hemelb::lb::collisions::Collision* LBM::GetCollision(int i)
    {
  switch(i)
  {
    case 0: return mMidFluidCollision;
    case 1: return mWallCollision;
    case 2: return mInletCollision;
    case 3: return mOutletCollision;
    case 4: return mInletWallCollision;
    case 5: return mOutletWallCollision;
  }
  return NULL;
    }

// The entire simulation time step takes place through this function
// when the convergence criterion is not applied. Communications
// automatically handle the streaming stage pertaining to neighbouring
// subdomains.
int LBM::lbmCycle (int perform_rt, Net *net)
{
  double density, vx, vy, vz, velocity;
  int collision_type, offset;
  int i, m;

  NeighProc *neigh_proc_p;


  for (m = 0; m < net->neigh_procs; m++)
    {
      neigh_proc_p = &net->neigh_proc[ m ];
#ifndef NOMPI
      net->err = MPI_Irecv (&f_old[ neigh_proc_p->f_head ],
			    neigh_proc_p->fs, MPI_DOUBLE,
			    neigh_proc_p->id, 10, MPI_COMM_WORLD,
			    &net->req[ 0 ][ m ]);
#endif
    }

  offset = net->my_inner_sites;

  for (collision_type = 0; collision_type < COLLISION_TYPES; collision_type++)
  {
    for (i = offset; i < offset + net->my_inter_collisions[collision_type]; i++)
    {
      (*lbmUpdateSiteData[perform_rt]) (false, omega, i, &density, &vx, &vy, &vz, &velocity, net,
      this, GetCollision(collision_type));
    }
    offset += net->my_inter_collisions[collision_type];
  }

  for (m = 0; m < net->neigh_procs; m++)
    {
      neigh_proc_p = &net->neigh_proc[ m ];
#ifndef NOMPI
      net->err = MPI_Isend (&f_new[ neigh_proc_p->f_head ],
			    neigh_proc_p->fs, MPI_DOUBLE,
			    neigh_proc_p->id, 10, MPI_COMM_WORLD,
			    &net->req[ 0 ][ net->neigh_procs + m ]);
#endif
    }

  offset = 0;

  for (collision_type = 0; collision_type < COLLISION_TYPES; collision_type++)
    {
      for (i = offset; i < offset + net->my_inner_collisions[ collision_type ]; i++)
	{
        (*lbmUpdateSiteData[perform_rt]) (false, omega, i, &density, &vx, &vy, &vz, &velocity, net,
        this, GetCollision(collision_type));
	}
      offset += net->my_inner_collisions[ collision_type ];
    }

  for (m = 0; m < net->neigh_procs; m++)
    {
#ifndef NOMPI
      net->err = MPI_Wait (&net->req[ 0 ][ m ], net->status);
      net->err = MPI_Wait (&net->req[ 0 ][ net->neigh_procs + m ], net->status);
#endif
    }

  // Copy the distribution functions received from the neighbouring
  // processors into the destination buffer "f_new".
  for (i = 0; i < net->shared_fs; i++)
    {
      f_new[ f_recv_iv[i] ] = f_old[ net->neigh_proc[0].f_head + i ];
    }

  // Do any cleanup steps necessary on boundary nodes
  // TODO: We should probably do this for all collision types
  // (Overhead on unnecessary calls will eventually be tiny, even if it isn't now.
  offset = net->my_inner_sites + net->my_inter_collisions[0];

  for (i = offset; i < offset + net->my_inter_collisions[ 1 ]; i++)
  {
    (*lbmUpdateSiteData[perform_rt]) (true, omega, i, &density, &vx, &vy, &vz, &velocity, net,
    this, GetCollision(1));
  }

  offset = net->my_inner_collisions[ 0 ];

  for (i = offset; i < offset + net->my_inner_collisions[ 1 ]; i++)
  {
    (*lbmUpdateSiteData[perform_rt]) (true, omega, i, &density, &vx, &vy, &vz, &velocity, net,
    this, GetCollision(1));
  }

  // Swap f_old and f_new ready for the next timestep.
  double *temp = f_old;
  f_old = f_new;
  f_new = temp;

  return STABLE;
}

void LBM::lbmCalculateFlowFieldValues ()
{
  double *local_data;
  double *global_data;

  int i;

  int lMaxInlets = hemelb::util::max(6+inlets,2*inlets);

  local_data = new double[lMaxInlets];
  global_data = new double[lMaxInlets];

#ifndef NOMPI
  local_data[0] = lbm_density_min;
  local_data[1] = lbm_velocity_min;
  local_data[2] = lbm_stress_min;
  local_data[3] = lbm_density_max;
  local_data[4] = lbm_velocity_max;
  local_data[5] = lbm_stress_max;

  memcpy (&local_data[6], lbm_peak_inlet_velocity, sizeof(double) * inlets);

  MPI_Reduce (&local_data[0], &global_data[0], 3, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce (&local_data[3], &global_data[3], 3+inlets, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  lbm_density_min  = global_data[0];
  lbm_velocity_min = global_data[1];
  lbm_stress_min   = global_data[2];
  lbm_density_max  = global_data[3];
  lbm_velocity_max = global_data[4];
  lbm_stress_max   = global_data[5];

  memcpy (lbm_peak_inlet_velocity, &global_data[6], sizeof(double) * inlets);

  for (i = 0; i < inlets; i++)
    {
      local_data[ i ] = lbm_average_inlet_velocity[ i ];
      local_data[ inlets+i ] = lbm_inlet_count[ i ];
    }
  MPI_Reduce (local_data, global_data, 2*inlets, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  for (i = 0; i < inlets; i++)
    {
      lbm_average_inlet_velocity[ i ] = global_data[ i ];
      lbm_inlet_count[ i ] = global_data[ inlets+i ];
    }
#endif // NOMPI

  delete [] global_data;
  delete [] local_data;

  for (i = 0; i < inlets; i++)
    {
      lbm_average_inlet_velocity[i] /= lbm_inlet_count[i];
      lbm_average_inlet_velocity[i] = lbmConvertVelocityToPhysicalUnits (lbm_average_inlet_velocity[i]);
      lbm_peak_inlet_velocity[i] = lbmConvertVelocityToPhysicalUnits (lbm_peak_inlet_velocity[i]);
    }

  period = period;

  inlets = inlets;
}


int LBM::IsUnstable(Net *net)
{
  int is_unstable, stability;

  is_unstable = 0;

  for (int i = 0; i < net->my_sites; i++)
  {
    for (int l = 0; l < D3Q15::NUMVECTORS; l++)
    {
      if (f_old[i * D3Q15::NUMVECTORS + l] < 0.)
      {
        is_unstable = 1;
      }
    }
  }

#ifndef NOMPI
  net->err = MPI_Allreduce(&is_unstable, &stability, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  is_unstable = stability;
#endif

  return is_unstable;
}


// Update peak and average inlet velocities local to the current subdomain. 
void LBM::lbmUpdateInletVelocities (int time_step, Net *net)
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
	  lbm_peak_inlet_velocity[ i ] = -1e+30;
	  lbm_average_inlet_velocity[ i ] = 0.;
	  lbm_inlet_count[ i ] = 0;
	}
    }

      c1 = 15;
      c2 = 0;

      offset = net->my_inner_collisions[ 0 ] + net->my_inner_collisions[ 1 ];

  for (i = offset; i < offset + net->my_inner_collisions[ 2 ]; i++)
    {
      D3Q15::CalculateDensityAndVelocity(&f_old[ i*c1+c2 ], &density, &vx, &vy, &vz);

      inlet_id = (net->net_site_data[i] & BOUNDARY_ID_MASK) >> BOUNDARY_ID_SHIFT;

      if (is_inlet_normal_available)
	{
	  vx *= lbm_inlet_normal[ 3*inlet_id+0 ];
	  vy *= lbm_inlet_normal[ 3*inlet_id+1 ];
	  vz *= lbm_inlet_normal[ 3*inlet_id+2 ];

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
      lbm_peak_inlet_velocity[ inlet_id ] = fmax(lbm_peak_inlet_velocity[ inlet_id ], velocity);
      lbm_average_inlet_velocity[ inlet_id ] += velocity;
      ++lbm_inlet_count[ inlet_id ];
    }
  offset = net->my_inner_sites + net->my_inter_collisions[ 0 ] + net->my_inter_collisions[ 1 ];

  for (i = offset; i < offset + net->my_inter_collisions[ 2 ]; i++)
    {
      D3Q15::CalculateDensityAndVelocity(&f_old[ i*c1+c2 ], &density, &vx, &vy, &vz);

      inlet_id = (net->net_site_data[i] & BOUNDARY_ID_MASK) >> BOUNDARY_ID_SHIFT;

      if (is_inlet_normal_available)
	{
	  vx *= lbm_inlet_normal[ 3*inlet_id+0 ];
	  vy *= lbm_inlet_normal[ 3*inlet_id+1 ];
	  vz *= lbm_inlet_normal[ 3*inlet_id+2 ];

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
      lbm_peak_inlet_velocity[ inlet_id ] = fmax(lbm_peak_inlet_velocity[ inlet_id ], velocity);
      lbm_average_inlet_velocity[ inlet_id ] += velocity;
      ++lbm_inlet_count[ inlet_id ];
    }
}

// In the case of instability, this function restart the simulation
// with twice as many time steps per period and update the parameters
// that depends on this change.
void LBM::lbmRestart (Net *net)
{
  int i;

  for (i = 0; i < inlets; i++)
    {
      inlet_density_avg[i] = lbmConvertPressureToPhysicalUnits (inlet_density_avg[i] * Cs2);
      inlet_density_amp[i] = lbmConvertPressureGradToPhysicalUnits (inlet_density_amp[i] * Cs2);
    }
  for (i = 0; i < outlets; i++)
    {
      outlet_density_avg[i] = lbmConvertPressureToPhysicalUnits (outlet_density_avg[i] * Cs2);
      outlet_density_amp[i] = lbmConvertPressureGradToPhysicalUnits (outlet_density_amp[i] * Cs2);
    }
  period *= 2;

  for (i = 0; i < inlets; i++)
    {
      inlet_density_avg[i] = lbmConvertPressureToLatticeUnits (inlet_density_avg[i]) / Cs2;
      inlet_density_amp[i] = lbmConvertPressureGradToLatticeUnits (inlet_density_amp[i]) / Cs2;
    }
  for (i = 0; i < outlets; i++)
    {
      outlet_density_avg[i] = lbmConvertPressureToLatticeUnits (outlet_density_avg[i]) / Cs2;
      outlet_density_amp[i] = lbmConvertPressureGradToLatticeUnits (outlet_density_amp[i]) / Cs2;
    }

  RecalculateTauViscosityOmega ();

  lbmSetInitialConditions (net);
}

double LBM::GetMinPhysicalPressure()
{
  return lbmConvertPressureToPhysicalUnits(lbm_density_min * Cs2);
}
double LBM::GetMaxPhysicalPressure()
{
  return lbmConvertPressureToPhysicalUnits(lbm_density_max * Cs2);
}
double LBM::GetMinPhysicalVelocity()
{
  return lbmConvertVelocityToPhysicalUnits(lbm_velocity_min);
}
double LBM::GetMaxPhysicalVelocity()
{
  return lbmConvertVelocityToPhysicalUnits(lbm_velocity_max);
}
double LBM::GetMinPhysicalStress()
{
  return lbmConvertStressToPhysicalUnits(lbm_stress_min);
}
double LBM::GetMaxPhysicalStress()
{
  return lbmConvertStressToPhysicalUnits(lbm_stress_max);
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
