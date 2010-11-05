/*! \file config.cc
 \brief In this file, the functions useful for the input/output are reported
 */
#include <limits.h>
#include <sstream>
#include <math.h>
#include <string.h>

#include "debug/Debugger.h"
#include "lb.h"
#include "net.h"
#include "utilityFunctions.h"
#include "io/XdrFileReader.h"
#include "io/AsciiFileWriter.h"

//using namespace std;

/*!
 this function reads the XDR configuration file but does not store the system
 and calculate some parameters
 */

void LBM::lbmReadConfig(Net *net,
                        hemelb::lb::GlobalLatticeData &bGlobalLatticeData)
{
  /* Read the config file written by the segtool.
   *
   * All values encoded using XDR format. Uses int, double and u_int.
   * 
   * System parameters:
   *   double stress_type
   *   int blocks_x
   *   int blocks_y
   *   int blocks_z
   *   int block_size 
   *
   * For each block (all blocks_x * blocks_y * blocks_z of them): 
   *
   *   int flag (indicates presence of non-solid sites in the block)
   *   
   *   If flag == 0 go to next block
   *
   *   Otherwise for each site in the block (all block_size^3):
   *
   *     u_int site_data -- this is a bit field which indicates site
   *     type (OR with SITE_TYPE_MASK to get bits zero and one; 00 =
   *     solid, 01 = fluid, 10 = inlet, 11 = outlet) or edgeness (set
   *     bit with PRESSURE_EDGE_MASK)
   * 
   *     If solid or simple fluid, go to next site
   *     
   *     If inlet or outlet (irrespective of edge state) {
   *       double boundary_normal[3]
   *       double boundary_dist
   *     }
   *
   *     If edge bit set {
   *       double wall_normal[3]
   *       double wall_dist
   *     }
   *     
   *     double mDistanceToWall[14]
   */

  FILE* xdrFile = fopen(mSimConfig->DataFilePath.c_str(), "r");

  char* lProcIdentifier = net->GetCurrentProcIdentifier();

  if (xdrFile == NULL)
  {
    fprintf(stderr, "Unable to open file %s [%s], exiting\n",
            mSimConfig->DataFilePath.c_str(), lProcIdentifier);
    fflush(0x0);
    exit(0x0);
  }
  else
  {
    fprintf(stderr, "Opened config file %s [%s]\n",
            mSimConfig->DataFilePath.c_str(), lProcIdentifier);
  }
  fflush(NULL);

  delete[] lProcIdentifier;

  hemelb::io::XdrReader myReader = hemelb::io::XdrFileReader(xdrFile);

  int j, k, ii, jj, kk, l, m, n;
  int flag;

  unsigned int site_i, site_j, site_k;
  unsigned int *site_type;

  myReader.readDouble(lbm_stress_type);
  myReader.readInt(bGlobalLatticeData.BlocksX);
  myReader.readInt(bGlobalLatticeData.BlocksY);
  myReader.readInt(bGlobalLatticeData.BlocksZ);
  myReader.readInt(bGlobalLatticeData.BlockSize);

  bGlobalLatticeData.SitesX = bGlobalLatticeData.BlocksX
      * bGlobalLatticeData.BlockSize;
  bGlobalLatticeData.SitesY = bGlobalLatticeData.BlocksY
      * bGlobalLatticeData.BlockSize;
  bGlobalLatticeData.SitesZ = bGlobalLatticeData.BlocksZ
      * bGlobalLatticeData.BlockSize;

  bGlobalLatticeData.SitesPerBlockVolumeUnit = bGlobalLatticeData.BlockSize
      * bGlobalLatticeData.BlockSize * bGlobalLatticeData.BlockSize;

  // shift = log_2(block_size)
  int i = bGlobalLatticeData.BlockSize;
  bGlobalLatticeData.Log2BlockSize = 0;
  while (i > 1)
  {
    i >>= 1;
    ++bGlobalLatticeData.Log2BlockSize;
  }

  bGlobalLatticeData.BlockCount = bGlobalLatticeData.BlocksX
      * bGlobalLatticeData.BlocksY * bGlobalLatticeData.BlocksZ;

  net->map_block = new DataBlock[bGlobalLatticeData.BlockCount];

  total_fluid_sites = 0;

  site_min_x = INT_MAX;
  site_min_y = INT_MAX;
  site_min_z = INT_MAX;
  site_max_x = INT_MIN;
  site_max_y = INT_MIN;
  site_max_z = INT_MIN;

  net->fr_time = hemelb::util::myClock();

  n = -1;

  for (i = 0; i < bGlobalLatticeData.BlocksX; i++)
  {
    for (j = 0; j < bGlobalLatticeData.BlocksY; j++)
    {
      for (k = 0; k < bGlobalLatticeData.BlocksZ; k++)
      {
        ++n;

        net->map_block[n].site_data = NULL;
        net->map_block[n].ProcessorRankForEachBlockSite = NULL;
        net->map_block[n].wall_data = NULL;

        myReader.readInt(flag);

        if (flag == 0)
          continue;
        // Block contains some non-solid sites

        net->map_block[n].site_data
            = new unsigned int[bGlobalLatticeData.SitesPerBlockVolumeUnit];
        net->map_block[n].ProcessorRankForEachBlockSite
            = new int[bGlobalLatticeData.SitesPerBlockVolumeUnit];

        m = -1;

        for (ii = 0; ii < bGlobalLatticeData.BlockSize; ii++)
        {
          site_i = (i << bGlobalLatticeData.Log2BlockSize) + ii;

          for (jj = 0; jj < bGlobalLatticeData.BlockSize; jj++)
          {
            site_j = (j << bGlobalLatticeData.Log2BlockSize) + jj;

            for (kk = 0; kk < bGlobalLatticeData.BlockSize; kk++)
            {
              site_k = (k << bGlobalLatticeData.Log2BlockSize) + kk;

              ++m;

              site_type = &net->map_block[n].site_data[m];
              myReader.readUnsignedInt(*site_type);

              if ( (*site_type & SITE_TYPE_MASK) == hemelb::lb::SOLID_TYPE)
              {
                net->map_block[n].ProcessorRankForEachBlockSite[m] = 1 << 30;
                continue;
              }
              net->map_block[n].ProcessorRankForEachBlockSite[m] = -1;

              ++total_fluid_sites;

              site_min_x = hemelb::util::min(site_min_x, site_i);
              site_min_y = hemelb::util::min(site_min_y, site_j);
              site_min_z = hemelb::util::min(site_min_z, site_k);
              site_max_x = hemelb::util::max(site_max_x, site_i);
              site_max_y = hemelb::util::max(site_max_y, site_j);
              site_max_z = hemelb::util::max(site_max_z, site_k);

              if (lbm_stress_type == SHEAR_STRESS
                  && net->GetCollisionType(*site_type) != FLUID)
              {
                // Neither solid nor simple fluid
                if (net->map_block[n].wall_data == NULL)
                {
                  net->map_block[n].wall_data
                      = new WallData[bGlobalLatticeData.SitesPerBlockVolumeUnit];
                }

                if (net->GetCollisionType(*site_type) & INLET
                    || net->GetCollisionType(*site_type) & OUTLET)
                {
                  // INLET or OUTLET or both
                  for (l = 0; l < 3; l++)
                    myReader.readDouble(
                                        net->map_block[n].wall_data[m].boundary_nor[l]);

                  myReader.readDouble(
                                      net->map_block[n].wall_data[m].boundary_dist);
                }

                if (net->GetCollisionType(*site_type) & EDGE)
                {
                  // EDGE bit set
                  for (l = 0; l < 3; l++)
                    myReader.readDouble(
                                        net->map_block[n].wall_data[m].wall_nor[l]);

                  myReader.readDouble(net->map_block[n].wall_data[m].wall_dist);
                }

                for (l = 0; l < 14; l++)
                  myReader.readDouble(
                                      net->map_block[n].wall_data[m].cut_dist[l]);
              }
            } // kk
          } // jj
        } // ii
      } // k
    } // j
  } // i

  fclose(xdrFile);

  net->fr_time = hemelb::util::myClock() - net->fr_time;
}

/*!
 through this function the processor 0 reads the LB parameters
 and then communicate them to the other processors
 */
void LBM::lbmReadParameters(Net *net)
{

  double par_to_send[10000];
  int nParamsRead = 0;

  if (net->IsCurrentProcTheIOProc())
  {
    inlets = mSimConfig->Inlets.size();
    allocateInlets(inlets);

    for (int n = 0; n < inlets; n++)
    {
      hemelb::SimConfig::InOutLet *lInlet = mSimConfig->Inlets[n];

      inlet_density_avg[n] = lbmConvertPressureToLatticeUnits(lInlet->PMean)
          / Cs2;
      inlet_density_amp[n] = lbmConvertPressureGradToLatticeUnits(lInlet->PAmp)
          / Cs2;
      inlet_density_phs[n] = lInlet->PPhase * DEG_TO_RAD;
    }

    outlets = mSimConfig->Outlets.size();
    allocateOutlets(outlets);

    for (int n = 0; n < outlets; n++)
    {
      hemelb::SimConfig::InOutLet *lOutlet = mSimConfig->Outlets[n];
      outlet_density_avg[n] = lbmConvertPressureToLatticeUnits(lOutlet->PMean)
          / Cs2;
      outlet_density_amp[n]
          = lbmConvertPressureGradToLatticeUnits(lOutlet->PAmp) / Cs2;
      outlet_density_phs[n] = lOutlet->PPhase * DEG_TO_RAD;
    }

    lbm_average_inlet_velocity = new double[inlets];
    lbm_peak_inlet_velocity = new double[inlets];
    lbm_inlet_normal = new double[3 * inlets];
    lbm_inlet_count = new long int[inlets];

    is_inlet_normal_available = 1;

    for (int ii = 0; ii < inlets; ii++)
    {
      lbm_inlet_normal[3 * ii] = mSimConfig->Inlets[ii]->Normal.x;
      lbm_inlet_normal[3 * ii + 1] = mSimConfig->Inlets[ii]->Normal.y;
      lbm_inlet_normal[3 * ii + 2] = mSimConfig->Inlets[ii]->Normal.z;
    }

    par_to_send[0] = 0.1 + (double) inlets;
    par_to_send[1] = 0.1 + (double) outlets;
    par_to_send[2] = 0.1 + (double) is_inlet_normal_available;
  }
#ifndef NOMPI
  net->err = MPI_Bcast(par_to_send, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
  if (!net->IsCurrentProcTheIOProc())
  {
    inlets = (int) par_to_send[0];
    outlets = (int) par_to_send[1];
    is_inlet_normal_available = (int) par_to_send[2];

    allocateInlets(inlets);
    allocateOutlets(outlets);

    lbm_average_inlet_velocity = new double[inlets];
    lbm_peak_inlet_velocity = new double[inlets];
    lbm_inlet_normal = new double[3 * inlets];
    lbm_inlet_count = new long int[inlets];
  }
  else
  {
    for (int n = 0; n < inlets; n++)
    {
      par_to_send[3 * n + 0] = inlet_density_avg[n];
      par_to_send[3 * n + 1] = inlet_density_amp[n];
      par_to_send[3 * n + 2] = inlet_density_phs[n];
    }
    for (int n = 0; n < outlets; n++)
    {
      par_to_send[3 * inlets + 3 * n + 0] = outlet_density_avg[n];
      par_to_send[3 * inlets + 3 * n + 1] = outlet_density_amp[n];
      par_to_send[3 * inlets + 3 * n + 2] = outlet_density_phs[n];
    }
    if (is_inlet_normal_available)
    {
      for (int n = 0; n < inlets; n++)
      {
        par_to_send[3 * (inlets + outlets) + 3 * n + 0] = lbm_inlet_normal[3
            * n + 0];
        par_to_send[3 * (inlets + outlets) + 3 * n + 1] = lbm_inlet_normal[3
            * n + 1];
        par_to_send[3 * (inlets + outlets) + 3 * n + 2] = lbm_inlet_normal[3
            * n + 2];
      }
    }
  }
#ifndef NOMPI
  net->err = MPI_Bcast(par_to_send, 3 * (inlets + outlets + inlets),
                       MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
  if (!net->IsCurrentProcTheIOProc())
  {
    for (int n = 0; n < inlets; n++)
    {
      inlet_density_avg[n] = par_to_send[3 * n + 0];
      inlet_density_amp[n] = par_to_send[3 * n + 1];
      inlet_density_phs[n] = par_to_send[3 * n + 2];
    }
    for (int n = 0; n < outlets; n++)
    {
      outlet_density_avg[n] = par_to_send[3 * inlets + 3 * n + 0];
      outlet_density_amp[n] = par_to_send[3 * inlets + 3 * n + 1];
      outlet_density_phs[n] = par_to_send[3 * inlets + 3 * n + 2];
    }
    if (is_inlet_normal_available)
    {
      for (int n = 0; n < inlets; n++)
      {
        lbm_inlet_normal[3 * n + 0] = par_to_send[3 * (inlets + outlets) + 3
            * n + 0];
        lbm_inlet_normal[3 * n + 1] = par_to_send[3 * (inlets + outlets) + 3
            * n + 1];
        lbm_inlet_normal[3 * n + 2] = par_to_send[3 * (inlets + outlets) + 3
            * n + 2];
      }
    }
  }
  lbmUpdateBoundaryDensities(0, 0);

  RecalculateTauViscosityOmega();
}

void LBM::allocateInlets(int nInlets)
{
  nInlets = hemelb::util::max(1, nInlets);
  inlet_density = new double[nInlets];
  inlet_density_avg = new double[nInlets];
  inlet_density_amp = new double[nInlets];
  inlet_density_phs = new double[nInlets];
}

void LBM::allocateOutlets(int nOutlets)
{
  nOutlets = hemelb::util::max(1, nOutlets);
  outlet_density = new double[nOutlets];
  outlet_density_avg = new double[nOutlets];
  outlet_density_amp = new double[nOutlets];
  outlet_density_phs = new double[nOutlets];
}

void LBM::lbmWriteConfig(int stability,
                         std::string output_file_name,
                         Net *net,
                         hemelb::lb::GlobalLatticeData &iGlobalLatticeData,
                         hemelb::lb::LocalLatticeData &iLocalLatticeData)
{
  /* This routine writes the flow field on file. The data are gathered
   to the root processor and written from there.  The format
   comprises:

   0- Flag for simulation stability, 0 or 1

   1- Voxel size in physical units (units of m)

   2- vertex coords of the minimum bounding box with minimum values
   (x, y and z values)

   3- vertex coords of the minimum bounding box with maximum values
   (x, y and z values)

   4- #voxels within the minimum bounding box along the x, y, z axes
   (3 values)

   5- total number of fluid voxels

   6-And then a list of the fluid voxels... For each fluid voxel:

   a- the (x, y, z) coordinates in lattice units (3 values)
   b- the pressure in physical units (mmHg)
   c- (x,y,z) components of the velocity field in physical units (3
   values, m/s)
   d- the von Mises stress in physical units (Pa) (the stored shear
   stress is equal to -1 if the fluid voxel is not at the wall)

   */
  hemelb::io::AsciiFileWriter *realSnap = NULL;

  float *local_flow_field, *gathered_flow_field;

  double density;
  double pressure;
  double vx, vy, vz;
  double stress;
  double f_eq[D3Q15::NUMVECTORS], f_neq[D3Q15::NUMVECTORS];
  float pressure_par, velocity_par, stress_par;

  int buffer_size;
  int fluid_sites_max;
  int communication_period, communication_iters;
  int comPeriodDelta, iters;
  int par;
  int shrinked_sites_x, shrinked_sites_y, shrinked_sites_z;

  short int *local_site_data, *gathered_site_data;

  unsigned int my_site_id;

  // parameters useful to convert pressure, velocity and stress from
  // lattice to physical units
  pressure_par = PULSATILE_PERIOD / (period * voxel_size * voxel_size);

  pressure_par = BLOOD_DENSITY / (mmHg_TO_PASCAL * pressure_par * pressure_par
      * voxel_size * voxel_size);

  velocity_par = 1.0 / (voxel_size * ( (tau - 0.5) / 3.) / (BLOOD_VISCOSITY
      / BLOOD_DENSITY));

  stress_par = ( (tau - 0.5) / 3.0) / (BLOOD_VISCOSITY / BLOOD_DENSITY);
  stress_par = BLOOD_DENSITY / (stress_par * stress_par * voxel_size
      * voxel_size);

  if (net->IsCurrentProcTheIOProc())
  {
    realSnap = new hemelb::io::AsciiFileWriter(output_file_name);
    //snap << stability << snap->eol;
    (*realSnap << stability) << realSnap->eol;
    //snap->write(stability); snap->writeRecordSeparator();
  }
  hemelb::io::Writer& snap = *realSnap;

  if (stability == UNSTABLE)
  {
    if (net->IsCurrentProcTheIOProc())
    {
      delete realSnap;
    }
    return;
  }

  if (net->IsCurrentProcTheIOProc())
  {
    shrinked_sites_x = 1 + site_max_x - site_min_x;
    shrinked_sites_y = 1 + site_max_y - site_min_y;
    shrinked_sites_z = 1 + site_max_z - site_min_z;

    snap << voxel_size << hemelb::io::Writer::eol;
    snap << site_min_x << site_min_y << site_min_z << hemelb::io::Writer::eol;
    snap << site_max_x << site_max_y << site_max_z << hemelb::io::Writer::eol;
    snap << shrinked_sites_x << shrinked_sites_y << shrinked_sites_z
        << hemelb::io::Writer::eol;
    snap << total_fluid_sites << hemelb::io::Writer::eol;
  }

  fluid_sites_max = 0;

  for (int n = 0; n < net->mProcessorCount; n++)
  {
    fluid_sites_max = hemelb::util::max(fluid_sites_max,
                                        net->mFluidSitesOnEachProcessor[n]);
  }

  // "buffer_size" is the size of the flow field buffer to send to the
  // root processor ("local_flow_field") and that to accommodate the
  // received ones from the non-root processors
  // ("gathered_flow_field").  If "buffer_size" is larger the
  // frequency with which data communication to the root processor is
  // performed becomes lower and viceversa
  buffer_size = hemelb::util::min(1000000, fluid_sites_max
      * net->mProcessorCount);

  communication_period = int(ceil(double(buffer_size) / net->mProcessorCount));

  communication_iters = hemelb::util::max(1, int(ceil(double(fluid_sites_max)
      / communication_period)));

  local_flow_field = new float[MACROSCOPIC_PARS * communication_period];
  gathered_flow_field = new float[MACROSCOPIC_PARS * communication_period
      * net->mProcessorCount];

  local_site_data = new short int[3 * communication_period];
  gathered_site_data = new short int[3 * communication_period
      * net->mProcessorCount];

  for (comPeriodDelta = 0; comPeriodDelta < communication_period; comPeriodDelta++)
  {
    local_site_data[comPeriodDelta * 3] = -1;
  }
  iters = 0;
  comPeriodDelta = 0;

  par = 0;

  /* The following loops scan over every single macrocell (block). If
   the block is non-empty, it scans the fluid sites within that block
   If the site is fluid, it calculates the flow field and then is
   converted to physical units and stored in a buffer to send to the
   root processor */

  int n = -1; // net->proc_block counter
  for (int i = 0; i < iGlobalLatticeData.SitesX; i
      += iGlobalLatticeData.BlockSize)
  {
    for (int j = 0; j < iGlobalLatticeData.SitesY; j
        += iGlobalLatticeData.BlockSize)
    {
      for (int k = 0; k < iGlobalLatticeData.SitesZ; k
          += iGlobalLatticeData.BlockSize)
      {

        ++n;

        if (net->map_block[n].ProcessorRankForEachBlockSite == NULL)
        {
          continue;
        }
        int m = -1;

        for (int site_i = i; site_i < i + iGlobalLatticeData.BlockSize; site_i++)
        {
          for (int site_j = j; site_j < j + iGlobalLatticeData.BlockSize; site_j++)
          {
            for (int site_k = k; site_k < k + iGlobalLatticeData.BlockSize; site_k++)
            {

              m++;
              if (!net->IsCurrentProcRank(
                                          net->map_block[n].ProcessorRankForEachBlockSite[m]))
              {
                continue;
              }

              my_site_id = net->map_block[n].site_data[m];

              /* No idea what this does */
              if (my_site_id & (1U << 31U))
                continue;

              if (my_site_id == 31702)
                hemelb::debug::Debugger::Get()->BreakHere();

              // TODO Utter filth. The cases where the whole site data is exactly equal
              // to "FLUID_TYPE" and where just the type-component of the whole site data
              // is equal to "FLUID_TYPE" are handled differently.
              if (iLocalLatticeData.mSiteData[my_site_id]
                  == hemelb::lb::FLUID_TYPE)
              {
                D3Q15::CalculateDensityVelocityFEq(
                                                   &iLocalLatticeData.FOld[ (my_site_id
                                                       * (par + 1) + par)
                                                       * D3Q15::NUMVECTORS],
                                                   density, vx, vy, vz, f_eq);

                for (unsigned int l = 0; l < D3Q15::NUMVECTORS; l++)
                {
                  f_neq[l] = iLocalLatticeData.FOld[ (my_site_id * (par + 1)
                      + par) * D3Q15::NUMVECTORS + l] - f_eq[l];
                }

              }
              else
              { // not FLUID_TYPE
                lbmCalculateBC(&iLocalLatticeData.FOld[ (my_site_id * (par + 1)
                    + par) * D3Q15::NUMVECTORS],
                               iLocalLatticeData.GetSiteType(my_site_id),
                               iLocalLatticeData.GetBoundaryId(my_site_id),
                               &density, &vx, &vy, &vz, f_neq);
              }

              if (my_site_id == 31702)
                hemelb::debug::Debugger::Get()->BreakHere();

              if (lbm_stress_type == SHEAR_STRESS)
              {
                if (iLocalLatticeData.GetNormalToWall(my_site_id)[0] >= 1.0e+30)
                {
                  stress = -1.0;
                }
                else
                {
                  D3Q15::CalculateShearStress(
                                              density,
                                              f_neq,
                                              &iLocalLatticeData.GetNormalToWall(
                                                                                 my_site_id)[0],
                                              stress, lbm_stress_par);
                }
              }
              else
              {
                D3Q15::CalculateVonMisesStress(f_neq, stress, lbm_stress_par);
              }

              vx /= density;
              vy /= density;
              vz /= density;

              // conversion from lattice to physical units
              pressure = REFERENCE_PRESSURE + ( (density - 1.0) * Cs2)
                  * pressure_par;

              vx *= velocity_par;
              vy *= velocity_par;
              vz *= velocity_par;

              stress *= stress_par;

              local_flow_field[MACROSCOPIC_PARS * comPeriodDelta + 0]
                  = float(pressure);
              local_flow_field[MACROSCOPIC_PARS * comPeriodDelta + 1]
                  = float(vx);
              local_flow_field[MACROSCOPIC_PARS * comPeriodDelta + 2]
                  = float(vy);
              local_flow_field[MACROSCOPIC_PARS * comPeriodDelta + 3]
                  = float(vz);
              local_flow_field[MACROSCOPIC_PARS * comPeriodDelta + 4]
                  = float(stress);

              local_site_data[3 * comPeriodDelta + 0] = site_i;
              local_site_data[3 * comPeriodDelta + 1] = site_j;
              local_site_data[3 * comPeriodDelta + 2] = site_k;

              if (++comPeriodDelta != communication_period)
                continue;

              comPeriodDelta = 0;
              ++iters;
#ifndef NOMPI
              net->err = MPI_Gather(local_flow_field, MACROSCOPIC_PARS
                  * communication_period, MPI_FLOAT, gathered_flow_field,
                                    MACROSCOPIC_PARS * communication_period,
                                    MPI_FLOAT, 0, MPI_COMM_WORLD);

              net->err = MPI_Gather(local_site_data, 3 * communication_period,
                                    MPI_SHORT, gathered_site_data, 3
                                        * communication_period, MPI_SHORT, 0,
                                    MPI_COMM_WORLD);
#endif
              if (net->IsCurrentProcTheIOProc())
              {

                for (int l = 0; l < net->mProcessorCount * communication_period; l++)
                {
                  if (gathered_site_data[l * 3 + 0] == -1)
                    continue;

                  gathered_site_data[l * 3 + 0] -= site_min_x;
                  gathered_site_data[l * 3 + 1] -= site_min_y;
                  gathered_site_data[l * 3 + 2] -= site_min_z;

                  snap << gathered_site_data[l * 3 + 0] << gathered_site_data[l
                      * 3 + 1] << gathered_site_data[l * 3 + 2];

                  for (int kk = 0; kk < MACROSCOPIC_PARS; kk++)
                  {
                    snap << gathered_flow_field[MACROSCOPIC_PARS * l + kk];
                  }
                  snap << hemelb::io::Writer::eol;
                }

              }

              for (int l = 0; l < communication_period; l++)
              {
                local_site_data[l * 3] = -1;
              }

            } // for site_k
          } // for site_j
        } // for site_i

      } // for k
    } // for j
  } //for i

  if (iters != communication_iters)
  {
    ++iters;

    // Weirdly initialized for
    for (; iters <= communication_iters; iters++)
    {
#ifndef NOMPI
      net->err = MPI_Gather(local_flow_field, MACROSCOPIC_PARS
          * communication_period, MPI_FLOAT, gathered_flow_field,
                            MACROSCOPIC_PARS * communication_period, MPI_FLOAT,
                            0, MPI_COMM_WORLD);

      net->err = MPI_Gather(local_site_data, 3 * communication_period,
                            MPI_SHORT, gathered_site_data, 3
                                * communication_period, MPI_SHORT, 0,
                            MPI_COMM_WORLD);
#endif

      if (net->IsCurrentProcTheIOProc())
      {
        for (int l = 0; l < net->mProcessorCount * communication_period; l++)
        {

          if (gathered_site_data[l * 3 + 0] == -1)
            continue;

          gathered_site_data[l * 3 + 0] -= site_min_x;
          gathered_site_data[l * 3 + 1] -= site_min_y;
          gathered_site_data[l * 3 + 2] -= site_min_z;

          snap << gathered_site_data[l * 3 + 0]
              << gathered_site_data[l * 3 + 1] << gathered_site_data[l * 3 + 2];

          for (int kk = 0; kk < MACROSCOPIC_PARS; kk++)
          {
            snap << gathered_flow_field[MACROSCOPIC_PARS * l + kk];
          }
          snap << hemelb::io::Writer::eol;
        }
      }

      for (int l = 0; l < communication_period; l++)
      {
        local_site_data[l * 3] = -1;
      }

    } // weird for

  }

  if (net->IsCurrentProcTheIOProc())
  {
    delete realSnap;
  }

  delete[] gathered_site_data;
  delete[] local_site_data;
  delete[] gathered_flow_field;
  delete[] local_flow_field;
}

void LBM::ReadVisParameters(Net *net)
{
  float lDensity_threshold_min, lDensity_threshold_minmax_inv,
      lVelocity_threshold_max_inv, lStress_threshold_max_inv;
  float par_to_send[9];
  float density_min, density_max, velocity_max, stress_max;

  int i;

  if (net->IsCurrentProcTheIOProc())
  {
    velocity_max = lbmConvertVelocityToLatticeUnits(mSimConfig->MaxVelocity);
    stress_max = lbmConvertStressToLatticeUnits(mSimConfig->MaxStress);

    par_to_send[0] = mSimConfig->VisCentre.x;
    par_to_send[1] = mSimConfig->VisCentre.y;
    par_to_send[2] = mSimConfig->VisCentre.z;
    par_to_send[3] = mSimConfig->VisLongitude;
    par_to_send[4] = mSimConfig->VisLatitude;
    par_to_send[5] = mSimConfig->VisZoom;
    par_to_send[6] = mSimConfig->VisBrightness;
    par_to_send[7] = velocity_max;
    par_to_send[8] = stress_max;
  }
#ifndef NOMPI
  net->err = MPI_Bcast(par_to_send, 9, MPI_FLOAT, 0, MPI_COMM_WORLD);
#endif

  mSimConfig->VisCentre.x = par_to_send[0];
  mSimConfig->VisCentre.y = par_to_send[1];
  mSimConfig->VisCentre.z = par_to_send[2];
  mSimConfig->VisLongitude = par_to_send[3];
  mSimConfig->VisLatitude = par_to_send[4];
  mSimConfig->VisZoom = par_to_send[5];
  mSimConfig->VisBrightness = par_to_send[6];
  velocity_max = par_to_send[7];
  stress_max = par_to_send[8];

  density_min = +1.0e+30F;
  density_max = -1.0e+30F;

  for (i = 0; i < inlets; i++)
  {
    density_min = fminf(density_min, inlet_density_avg[i]
        - inlet_density_amp[i]);
    density_max = fmaxf(density_max, inlet_density_avg[i]
        + inlet_density_amp[i]);
  }
  for (i = 0; i < outlets; i++)
  {
    density_min = fminf(density_min, outlet_density_avg[i]
        - outlet_density_amp[i]);
    density_max = fmaxf(density_max, outlet_density_avg[i]
        + outlet_density_amp[i]);
  }
  lDensity_threshold_min = density_min;

  lDensity_threshold_minmax_inv = 1.0F / (density_max - density_min);
  lVelocity_threshold_max_inv = 1.0F / velocity_max;
  lStress_threshold_max_inv = 1.0F / stress_max;

  hemelb::vis::controller->SetSomeParams(mSimConfig->VisBrightness,
                                         lDensity_threshold_min,
                                         lDensity_threshold_minmax_inv,
                                         lVelocity_threshold_max_inv,
                                         lStress_threshold_max_inv);
}

