/*! \file config.cc
 \brief In this file, the functions useful for the input/output are reported
 */
#include <limits.h>
#include <sstream>
#include <math.h>
#include <string.h>

#include "debug/Debugger.h"
#include "lb/lb.h"
#include "net/net.h"
#include "util/utilityFunctions.h"
#include "io/XdrMemReader.h"
#include "io/XdrMemWriter.h"
#include "io/AsciiFileWriter.h"
#include "geometry/LatticeData.h"

namespace hemelb
{
  namespace lb
  {
    void LBM::ReadParameters()
    {
      inlets = (int) mSimConfig->Inlets.size();
      outlets = (int) mSimConfig->Outlets.size();

      inlet_normal = new distribn_t[3 * inlets];

      for (int ii = 0; ii < inlets; ii++)
      {
        inlet_normal[3 * ii] = mSimConfig->Inlets[ii].Normal.x;
        inlet_normal[3 * ii + 1] = mSimConfig->Inlets[ii].Normal.y;
        inlet_normal[3 * ii + 2] = mSimConfig->Inlets[ii].Normal.z;
      }

      RecalculateTauViscosityOmega();
    }

    void LBM::WriteConfigParallel(hemelb::lb::Stability stability,
                                  std::string output_file_name,
                                  BoundaryComms* iInletComms,
                                  BoundaryComms* iOutletComms)
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

      if (stability == hemelb::lb::Unstable)
      {
        MPI_File_delete(&output_file_name[0], MPI_INFO_NULL);
        return;
      }

      MPI_Status lStatus;

      MPI_File lOutputFile;

      MPI_File_open(MPI_COMM_WORLD,
                    &output_file_name[0],
                    MPI_MODE_WRONLY | MPI_MODE_CREATE,
                    MPI_INFO_NULL,
                    &lOutputFile);

      /* Preamble has an enum (int) for stability, a double for voxel size,
       * 3 ints for minimum (x,y,z) in bounding box, 3 ints for maximum (x,y,z)
       * in bounding box, 3 ints for number of coords in each of (x,y,z),
       * 1 int for number of fluid voxels.*/
      const int lPreambleLength = 4 + 8 + (3 * 4) + (3 * 4) + (3 * 4) + 4;

      std::string lReadMode = "native";

      MPI_Datatype viewType = MpiDataType<char> ();
      MPI_File_set_view(lOutputFile, 0, viewType, viewType, &lReadMode[0], MPI_INFO_NULL);

      topology::NetworkTopology* netTop = topology::NetworkTopology::Instance();

      if (netTop->IsCurrentProcTheIOProc())
      {
        char lBuffer[lPreambleLength];
        hemelb::io::XdrMemWriter lWriter = hemelb::io::XdrMemWriter(lBuffer, lPreambleLength);

        lWriter << stability << voxel_size << (int) siteMins[0] << (int) siteMins[1]
            << (int) siteMins[2] << (int) siteMaxes[0] << (int) siteMaxes[1] << (int) siteMaxes[2]
            << (int) (1 + siteMaxes[0] - siteMins[0]) << (int) (1 + siteMaxes[1] - siteMins[1])
            << (int) (1 + siteMaxes[2] - siteMins[2]) << (int) total_fluid_sites;

        MPI_File_write(lOutputFile, lBuffer, lPreambleLength, MpiDataType(lBuffer[0]), &lStatus);
      }

      /*
       For each fluid voxel, we write
       a- the (x, y, z) coordinates in lattice units (3 ints)
       b- the pressure in physical units (mmHg, 1 x float)
       c- (x,y,z) components of the velocity field in physical units (3
       values, m/s, floats)
       d- the von Mises stress in physical units (Pa) (the stored shear
       stress is equal to -1 if the fluid voxel is not at the wall, 1 x float)
       */

      const short int lOneFluidSiteLength = (3 * 4) + (5 * 4);

      site_t lLocalSitesInitialOffset = lPreambleLength;

      for (proc_t ii = 0; ii < netTop->GetLocalRank(); ii++)
      {
        lLocalSitesInitialOffset += lOneFluidSiteLength * netTop->FluidSitesOnEachProcessor[ii];
      }

      MPI_File_set_view(lOutputFile,
                        lLocalSitesInitialOffset,
                        viewType,
                        viewType,
                        &lReadMode[0],
                        MPI_INFO_NULL);

      site_t lLocalWriteLength = lOneFluidSiteLength
          * netTop->FluidSitesOnEachProcessor[netTop->GetLocalRank()];
      char * lFluidSiteBuffer = new char[lLocalWriteLength];
      hemelb::io::XdrMemWriter lWriter = hemelb::io::XdrMemWriter(lFluidSiteBuffer,
                                                                  (unsigned int) lLocalWriteLength);

      /* The following loops scan over every single macrocell (block). If
       the block is non-empty, it scans the fluid sites within that block
       If the site is fluid, it calculates the flow field and then is
       converted to physical units and stored in a buffer to send to the
       root processor */

      site_t n = -1;
      for (site_t i = 0; i < mLatDat->GetXSiteCount(); i += mLatDat->GetBlockSize())
      {
        for (site_t j = 0; j < mLatDat->GetYSiteCount(); j += mLatDat->GetBlockSize())
        {
          for (site_t k = 0; k < mLatDat->GetZSiteCount(); k += mLatDat->GetBlockSize())
          {

            ++n;

            if (mLatDat->GetBlock(n)->ProcessorRankForEachBlockSite == NULL)
            {
              continue;
            }
            site_t m = -1;

            for (site_t site_i = i; site_i < i + mLatDat->GetBlockSize(); site_i++)
            {
              for (site_t site_j = j; site_j < j + mLatDat->GetBlockSize(); site_j++)
              {
                for (site_t site_k = k; site_k < k + mLatDat->GetBlockSize(); site_k++)
                {

                  m++;
                  if (netTop->GetLocalRank()
                      != mLatDat->GetBlock(n)->ProcessorRankForEachBlockSite[m])
                  {
                    continue;
                  }

                  site_t my_site_id = mLatDat->GetBlock(n)->site_data[m];

                  /* No idea what this does */
                  if (my_site_id & BIG_NUMBER3)
                    continue;

                  distribn_t density, vx, vy, vz, f_eq[D3Q15::NUMVECTORS],
                      f_neq[D3Q15::NUMVECTORS], stress, pressure;

                  // TODO Utter filth. The cases where the whole site data is exactly equal
                  // to "FLUID_TYPE" and where just the type-component of the whole site data
                  // is equal to "FLUID_TYPE" are handled differently.
                  if (mLatDat->GetSiteData(my_site_id) == geometry::LatticeData::FLUID_TYPE)
                  {
                    D3Q15::CalculateDensityVelocityFEq(mLatDat->GetFOld(my_site_id
                        * D3Q15::NUMVECTORS), density, vx, vy, vz, f_eq);

                    for (unsigned int l = 0; l < D3Q15::NUMVECTORS; l++)
                    {
                      f_neq[l] = *mLatDat->GetFOld(my_site_id * D3Q15::NUMVECTORS + l) - f_eq[l];
                    }

                  }
                  else
                  { // not FLUID_TYPE
                    mBoundaryValues->CalculateBC(mLatDat->GetFOld(my_site_id * D3Q15::NUMVECTORS),
                                                 mLatDat->GetSiteType(my_site_id),
                                                 mLatDat->GetBoundaryId(my_site_id),
                                                 &density,
                                                 &vx,
                                                 &vy,
                                                 &vz,
                                                 f_neq,
                                                 iInletComms,
                                                 iOutletComms);
                  }

                  if (mParams.StressType == hemelb::lb::ShearStress)
                  {
                    if (mLatDat->GetNormalToWall(my_site_id)[0] >= NO_VALUE)
                    {
                      stress = -1.0;
                    }
                    else
                    {
                      D3Q15::CalculateShearStress(density,
                                                  f_neq,
                                                  &mLatDat->GetNormalToWall(my_site_id)[0],
                                                  stress,
                                                  mParams.StressParameter);
                    }
                  }
                  else
                  {
                    D3Q15::CalculateVonMisesStress(f_neq, stress, mParams.StressParameter);
                  }

                  vx /= density;
                  vy /= density;
                  vz /= density;

                  // conversion from lattice to physical units
                  pressure = mUnits->ConvertPressureToPhysicalUnits(density * Cs2);

                  vx = mUnits->ConvertVelocityToPhysicalUnits(vx);
                  vy = mUnits->ConvertVelocityToPhysicalUnits(vy);
                  vz = mUnits->ConvertVelocityToPhysicalUnits(vz);

                  stress = mUnits->ConvertStressToPhysicalUnits(stress);

                  lWriter << (int) (site_i - siteMins[0]) << (int) (site_j - siteMins[1])
                      << (int) (site_k - siteMins[2]);

                  lWriter << float (pressure) << float (vx) << float (vy) << float (vz)
                      << float (stress);
                }
              }
            }
          }
        }
      }

      MPI_File_write_all(lOutputFile,
                         lFluidSiteBuffer,
                         (int) lLocalWriteLength,
                         MpiDataType(lFluidSiteBuffer[0]),
                         &lStatus);

      MPI_File_close(&lOutputFile);

      delete[] lFluidSiteBuffer;
    }

    void LBM::ReadVisParameters()
    {
      distribn_t density_min = std::numeric_limits<distribn_t>::max();
      distribn_t density_max = std::numeric_limits<distribn_t>::min();

      distribn_t velocity_max = mUnits->ConvertVelocityToLatticeUnits(mSimConfig->MaxVelocity);
      distribn_t stress_max = mUnits->ConvertStressToLatticeUnits(mSimConfig->MaxStress);

      for (int i = 0; i < inlets; i++)
      {
        density_min = util::NumericalFunctions::min(density_min,
                                                    mBoundaryValues->GetInletDensityMin(i));
        density_max = util::NumericalFunctions::max(density_max,
                                                    mBoundaryValues->GetInletDensityMax(i));
      }
      for (int i = 0; i < outlets; i++)
      {
        density_min = util::NumericalFunctions::min(density_min,
                                                    mBoundaryValues->GetOutletDensityMin(i));
        density_max = util::NumericalFunctions::max(density_max,
                                                    mBoundaryValues->GetOutletDensityMax(i));
      }

      distribn_t lDensity_threshold_min = density_min;
      distribn_t lDensity_threshold_minmax_inv = 1.0F / (density_max - density_min);
      distribn_t lVelocity_threshold_max_inv = 1.0F / velocity_max;
      distribn_t lStress_threshold_max_inv = 1.0F / stress_max;

      mVisControl->SetSomeParams(mSimConfig->VisBrightness,
                                 lDensity_threshold_min,
                                 lDensity_threshold_minmax_inv,
                                 lVelocity_threshold_max_inv,
                                 lStress_threshold_max_inv);
    }
  }
}
