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
#include "io/writers/xdr/XdrMemReader.h"
#include "io/writers/xdr/XdrMemWriter.h"
#include "io/writers/ascii/AsciiFileWriter.h"
#include "geometry/LatticeData.h"
#include "io/formats/snapshot.h"

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
        inlet_normal[3 * ii] = mSimConfig->Inlets[ii]->Normal.x;
        inlet_normal[3 * ii + 1] = mSimConfig->Inlets[ii]->Normal.y;
        inlet_normal[3 * ii + 2] = mSimConfig->Inlets[ii]->Normal.z;
      }
    }

    void LBM::WriteConfigParallel(hemelb::lb::Stability const stability,
                                  std::string output_file_name) const
    {
      /* This routine writes the flow field on file, using MPIO to coordinate
       * the writing. The format is detailed in io/formats/snapshot.h
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

      std::string lReadMode = "native";

      MPI_Datatype viewType = MpiDataType<char> ();
      MPI_File_set_view(lOutputFile, 0, viewType, viewType, &lReadMode[0], MPI_INFO_NULL);

      topology::NetworkTopology* netTop = topology::NetworkTopology::Instance();

      if (netTop->IsCurrentProcTheIOProc())
      {
        // Write the header according to format detailed in snapshot.h
        char lBuffer[io::formats::snapshot::HeaderLength];
        io::writers::xdr::XdrMemWriter lWriter =
            io::writers::xdr::XdrMemWriter(lBuffer, io::formats::snapshot::HeaderLength);
        lWriter << (unsigned int) io::formats::HemeLbMagicNumber
            << (unsigned int) io::formats::snapshot::MagicNumber
            << (unsigned int) io::formats::snapshot::VersionNumber;
        lWriter << (unsigned int) io::formats::snapshot::HeaderLength;
        lWriter << stability;
        lWriter << mLatDat->GetVoxelSize();
        lWriter << mLatDat->GetXOrigin() << mLatDat->GetYOrigin() << mLatDat->GetZOrigin();
        lWriter << (int) mLatDat->GetGlobalSiteMins().x << (int) mLatDat->GetGlobalSiteMins().y
            << (int) mLatDat->GetGlobalSiteMins().z;
        lWriter << (int) mLatDat->GetGlobalSiteMaxes().x << (int) mLatDat->GetGlobalSiteMaxes().y
            << (int) mLatDat->GetGlobalSiteMaxes().z;
        lWriter << (int) mLatDat->GetTotalFluidSites();

        MPI_File_write(lOutputFile,
                       lBuffer,
                       io::formats::snapshot::HeaderLength,
                       MpiDataType(lBuffer[0]),
                       &lStatus);
      }

      /* Now we write a record for each voxel.
       * Each task is responsible for creating (locally) a buffer which
       * contains all the records for the fluid sites for which it's
       * responsible. We then use MPIO to write these buffers (in rank order)
       * into the file after the header.
       */

      // This is the position in the file where the local task's voxels will be written.
      site_t lLocalSitesInitialOffset = io::formats::snapshot::HeaderLength;

      for (proc_t ii = 0; ii < netTop->GetLocalRank(); ii++)
      {
        lLocalSitesInitialOffset += io::formats::snapshot::VoxelRecordLength
            * mLatDat->GetFluidSiteCountOnProc(ii);
      }

      MPI_File_set_view(lOutputFile,
                        lLocalSitesInitialOffset,
                        viewType,
                        viewType,
                        &lReadMode[0],
                        MPI_INFO_NULL);

      site_t lLocalWriteLength = io::formats::snapshot::VoxelRecordLength
          * mLatDat->GetFluidSiteCountOnProc(netTop->GetLocalRank());
      char * lFluidSiteBuffer = new char[lLocalWriteLength];
      hemelb::io::writers::xdr::XdrMemWriter
          lWriter = hemelb::io::writers::xdr::XdrMemWriter(lFluidSiteBuffer,
                                                           (unsigned int) lLocalWriteLength);

      /* The following loops scan over every single macrocell (block). If
       * the block is non-empty, it scans the sites within that block. If the
       * site is fluid and present on the current task, it calculates the
       * flow field and encodes it to the local buffer.
       */

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
                                                           * D3Q15::NUMVECTORS),
                                                       density,
                                                       vx,
                                                       vy,
                                                       vz,
                                                       f_eq);

                    for (unsigned int l = 0; l < D3Q15::NUMVECTORS; l++)
                    {
                      f_neq[l] = *mLatDat->GetFOld(my_site_id * D3Q15::NUMVECTORS + l) - f_eq[l];
                    }

                  }
                  else
                  { // not FLUID_TYPE
                    CalculateBC(mLatDat->GetFOld(my_site_id * D3Q15::NUMVECTORS),
                                mLatDat->GetSiteType(my_site_id),
                                mLatDat->GetBoundaryId(my_site_id),
                                &density,
                                &vx,
                                &vy,
                                &vz,
                                f_neq);
                  }

                  if (mParams.StressType == hemelb::lb::ShearStress)
                  {
                    if ((mLatDat->GetSiteData(my_site_id) & PRESSURE_EDGE_MASK) == 0)
                    {
                      stress = -1.0;
                    }
                    else
                    {
                      D3Q15::CalculateShearStress(density,
                                                  f_neq,
                                                  &mLatDat->GetNormalToWall(my_site_id)[0],
                                                  stress,
                                                  mParams.GetStressParameter());
                    }
                  }
                  else
                  {
                    D3Q15::CalculateVonMisesStress(f_neq, stress, mParams.GetStressParameter());
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

                  const util::Vector3D<site_t>& siteMins = mLatDat->GetGlobalSiteMins();

                  lWriter << (int) (site_i - siteMins.x) << (int) (site_j - siteMins.y)
                      << (int) (site_k - siteMins.z);

                  lWriter << float(pressure) << float(vx) << float(vy) << float(vz)
                      << float(stress);
                }
              }
            }
          }
        }
      }
      // Hand the buffers over to MPIO to write to the file.
      MPI_File_write_all(lOutputFile,
                         lFluidSiteBuffer,
                         (int) lLocalWriteLength,
                         MpiDataType(lFluidSiteBuffer[0]),
                         &lStatus);

      MPI_File_close(&lOutputFile);

      delete[] lFluidSiteBuffer;
    }

    // Calculate the BCs for each boundary site type and the
    // non-equilibrium distribution functions.
    void LBM::CalculateBC(distribn_t f[],
                          hemelb::geometry::LatticeData::SiteType const iSiteType,
                          unsigned int const iBoundaryId,
                          distribn_t *density,
                          distribn_t *vx,
                          distribn_t *vy,
                          distribn_t *vz,
                          distribn_t f_neq[]) const
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
          *density = mInletValues->GetBoundaryDensity(iBoundaryId);
        }
        else
        {
          *density = mOutletValues->GetBoundaryDensity(iBoundaryId);
        }

        D3Q15::CalculateDensityAndVelocity(f, dummy_density, *vx, *vy, *vz);
        D3Q15::CalculateFeq(*density, *vx, *vy, *vz, f);

      }
      for (unsigned int l = 0; l < D3Q15::NUMVECTORS; l++)
      {
        f_neq[l] -= f[l];
      }

    }

    void LBM::ReadVisParameters()
    {
      distribn_t density_min = std::numeric_limits<distribn_t>::max();
      distribn_t density_max = std::numeric_limits<distribn_t>::min();

      distribn_t velocity_max = mUnits->ConvertVelocityToLatticeUnits(mSimConfig->MaxVelocity);
      distribn_t stress_max = mUnits->ConvertStressToLatticeUnits(mSimConfig->MaxStress);

      for (int i = 0; i < inlets; i++)
      {
        density_min = util::NumericalFunctions::min(density_min, mInletValues->GetDensityMin(i));
        density_max = util::NumericalFunctions::max(density_max, mInletValues->GetDensityMax(i));
      }
      for (int i = 0; i < outlets; i++)
      {
        density_min = util::NumericalFunctions::min(density_min, mOutletValues->GetDensityMin(i));
        density_max = util::NumericalFunctions::max(density_max, mOutletValues->GetDensityMax(i));
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
