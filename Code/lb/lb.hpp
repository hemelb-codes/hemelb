#ifndef HEMELB_LB_LB_HPP
#define HEMELB_LB_LB_HPP

#include "io/formats/snapshot.h"
#include "io/writers/xdr/XdrMemWriter.h"
#include "lb/lb.h"

namespace hemelb
{
  namespace lb
  {

    template<class LatticeType>
    hemelb::lb::LbmParameters* LBM<LatticeType>::GetLbmParams()
    {
      return &mParams;
    }

    template<class LatticeType>
    LBM<LatticeType>::LBM(configuration::SimConfig *iSimulationConfig,
                          net::Net* net,
                          geometry::LatticeData* latDat,
                          SimulationState* simState,
                          reporting::Timer &atimer) :
        mSimConfig(iSimulationConfig), mNet(net), mLatDat(latDat), mState(simState), mParams(PULSATILE_PERIOD_s
                                                                                                 / (distribn_t) simState->GetTimeStepsPerCycle(),
                                                                                             latDat->GetVoxelSize()), timer(atimer)
    {
      ReadParameters();
    }

    template<class LatticeType>
    void LBM<LatticeType>::CalculateMouseFlowField(const ScreenDensity densityIn,
                                                   const ScreenStress stressIn,
                                                   const LatticeDensity density_threshold_min,
                                                   const LatticeDensity density_threshold_minmax_inv,
                                                   const LatticeStress stress_threshold_max_inv,
                                                   PhysicalPressure &mouse_pressure,
                                                   PhysicalStress &mouse_stress)
    {
      LatticeDensity density = density_threshold_min + densityIn / density_threshold_minmax_inv;
      LatticeStress stress = stressIn / stress_threshold_max_inv;

      mouse_pressure = mUnits->ConvertPressureToPhysicalUnits(density * Cs2);
      mouse_stress = mUnits->ConvertStressToPhysicalUnits(stress);
    }

    template<class LatticeType>
    void LBM<LatticeType>::InitCollisions()
    {
      // TODO Note that the convergence checking is not yet implemented in the
      // new boundary condition hierarchy system.
      // It'd be nice to do this with something like
      // MidFluidCollision = new ConvergenceCheckingWrapper(new WhateverMidFluidCollision());

      kernels::InitParams initParams = kernels::InitParams();
      initParams.latDat = mLatDat;
      initParams.lbmParams = &mParams;

      initParams.siteCount = mLatDat->GetInnerCollisionCount(0) + mLatDat->GetInterCollisionCount(0);
      mMidFluidCollision = new tMidFluidCollision(initParams);

      initParams.siteCount = mLatDat->GetInnerCollisionCount(1) + mLatDat->GetInterCollisionCount(1);
      mWallCollision = new tWallCollision(initParams);

      initParams.siteCount = mLatDat->GetInnerCollisionCount(2) + mLatDat->GetInterCollisionCount(2);
      initParams.boundaryObject = mInletValues;
      mInletCollision = new tInletOutletCollision(initParams);

      initParams.siteCount = mLatDat->GetInnerCollisionCount(3) + mLatDat->GetInterCollisionCount(3);
      initParams.boundaryObject = mOutletValues;
      mOutletCollision = new tInletOutletCollision(initParams);

      initParams.siteCount = mLatDat->GetInnerCollisionCount(4) + mLatDat->GetInterCollisionCount(4);
      initParams.boundaryObject = mInletValues;
      mInletWallCollision = new tInletOutletWallCollision(initParams);

      initParams.siteCount = mLatDat->GetInnerCollisionCount(5) + mLatDat->GetInterCollisionCount(5);
      initParams.boundaryObject = mOutletValues;
      mOutletWallCollision = new tInletOutletWallCollision(initParams);
    }

    template<class LatticeType>
    void LBM<LatticeType>::Initialise(vis::Control* iControl,
                                      boundaries::BoundaryValues* iInletValues,
                                      boundaries::BoundaryValues* iOutletValues,
                                      util::UnitConverter* iUnits)
    {
      mInletValues = iInletValues;
      mOutletValues = iOutletValues;
      mUnits = iUnits;

      InitCollisions();

      SetInitialConditions();

      mVisControl = iControl;
    }

    template<class LatticeType>
    void LBM<LatticeType>::SetInitialConditions()
    {
      distribn_t density = 0.0;

      for (int i = 0; i < outlets; i++)
      {
        density += mOutletValues->GetDensityMin(i);
      }

      density /= outlets;

      for (site_t i = 0; i < mLatDat->GetLocalFluidSiteCount(); i++)
      {
        distribn_t f_eq[LatticeType::NUMVECTORS];

        LatticeType::CalculateFeq(density, 0.0, 0.0, 0.0, f_eq);

        geometry::Site site = mLatDat->GetSite(i);

        distribn_t* f_old_p = site.GetFOld();
        distribn_t* f_new_p = mLatDat->GetFNew(i * LatticeType::NUMVECTORS);

        for (unsigned int l = 0; l < LatticeType::NUMVECTORS; l++)
        {
          f_new_p[l] = f_old_p[l] = f_eq[l];
        }
      }
    }

    template<class LatticeType>
    void LBM<LatticeType>::RequestComms()
    {
      timer.Start();

      mLatDat->SendAndReceive(mNet);

      timer.Stop();
    }

    template<class LatticeType>
    void LBM<LatticeType>::PreSend()
    {
      timer.Start();

      site_t offset = mLatDat->GetInnerSiteCount();

      StreamAndCollide(mMidFluidCollision, offset, mLatDat->GetInterCollisionCount(0));
      offset += mLatDat->GetInterCollisionCount(0);

      StreamAndCollide(mWallCollision, offset, mLatDat->GetInterCollisionCount(1));
      offset += mLatDat->GetInterCollisionCount(1);

      mInletValues->FinishReceive();
      StreamAndCollide(mInletCollision, offset, mLatDat->GetInterCollisionCount(2));
      offset += mLatDat->GetInterCollisionCount(2);

      mOutletValues->FinishReceive();
      StreamAndCollide(mOutletCollision, offset, mLatDat->GetInterCollisionCount(3));
      offset += mLatDat->GetInterCollisionCount(3);

      StreamAndCollide(mInletWallCollision, offset, mLatDat->GetInterCollisionCount(4));
      offset += mLatDat->GetInterCollisionCount(4);

      StreamAndCollide(mOutletWallCollision, offset, mLatDat->GetInterCollisionCount(5));

      timer.Stop();
    }

    template<class LatticeType>
    void LBM<LatticeType>::PreReceive()
    {
      timer.Start();

      site_t offset = 0;

      StreamAndCollide(mMidFluidCollision, offset, mLatDat->GetInnerCollisionCount(0));
      offset += mLatDat->GetInnerCollisionCount(0);

      StreamAndCollide(mWallCollision, offset, mLatDat->GetInnerCollisionCount(1));
      offset += mLatDat->GetInnerCollisionCount(1);

      StreamAndCollide(mInletCollision, offset, mLatDat->GetInnerCollisionCount(2));
      offset += mLatDat->GetInnerCollisionCount(2);

      StreamAndCollide(mOutletCollision, offset, mLatDat->GetInnerCollisionCount(3));
      offset += mLatDat->GetInnerCollisionCount(3);

      StreamAndCollide(mInletWallCollision, offset, mLatDat->GetInnerCollisionCount(4));
      offset += mLatDat->GetInnerCollisionCount(4);

      StreamAndCollide(mOutletWallCollision, offset, mLatDat->GetInnerCollisionCount(5));

      timer.Stop();
    }

    template<class LatticeType>
    void LBM<LatticeType>::PostReceive()
    {
      timer.Start();

      // Copy the distribution functions received from the neighbouring
      // processors into the destination buffer "f_new".
      mLatDat->CopyReceived();

      // Do any cleanup steps necessary on boundary nodes
      site_t offset = 0;

      //TODO yup, this is horrible. If you read this, please improve the following code.
      PostStep(mMidFluidCollision, offset, mLatDat->GetInterCollisionCount(0));
      offset += mLatDat->GetInterCollisionCount(0);

      PostStep(mWallCollision, offset, mLatDat->GetInterCollisionCount(1));
      offset += mLatDat->GetInterCollisionCount(1);

      PostStep(mInletCollision, offset, mLatDat->GetInterCollisionCount(2));
      offset += mLatDat->GetInterCollisionCount(2);

      PostStep(mOutletCollision, offset, mLatDat->GetInterCollisionCount(3));
      offset += mLatDat->GetInterCollisionCount(3);

      PostStep(mInletWallCollision, offset, mLatDat->GetInterCollisionCount(4));
      offset += mLatDat->GetInterCollisionCount(4);

      PostStep(mOutletWallCollision, offset, mLatDat->GetInterCollisionCount(5));
      offset += mLatDat->GetInterCollisionCount(5);

      PostStep(mMidFluidCollision, offset, mLatDat->GetInnerCollisionCount(0));
      offset += mLatDat->GetInnerCollisionCount(0);

      PostStep(mWallCollision, offset, mLatDat->GetInnerCollisionCount(1));
      offset += mLatDat->GetInnerCollisionCount(1);

      PostStep(mInletCollision, offset, mLatDat->GetInnerCollisionCount(2));
      offset += mLatDat->GetInnerCollisionCount(2);

      PostStep(mOutletCollision, offset, mLatDat->GetInnerCollisionCount(3));
      offset += mLatDat->GetInnerCollisionCount(3);

      PostStep(mInletWallCollision, offset, mLatDat->GetInnerCollisionCount(4));
      offset += mLatDat->GetInnerCollisionCount(4);

      PostStep(mOutletWallCollision, offset, mLatDat->GetInnerCollisionCount(5));

      timer.Stop();
    }

    template<class LatticeType>
    void LBM<LatticeType>::EndIteration()
    {
      timer.Start();

      // Swap f_old and f_new ready for the next timestep.
      mLatDat->SwapOldAndNew();

      timer.Stop();
    }

    // In the case of instability, this function restart the simulation
    // with twice as many time steps per period and update the parameters
    // that depends on this change.
    template<class LatticeType>
    void LBM<LatticeType>::Reset()
    {
      mState->DoubleTimeResolution();

      mParams.Update(PULSATILE_PERIOD_s / (distribn_t) mState->GetTimeStepsPerCycle(), mLatDat->GetVoxelSize());

      SetInitialConditions();

      kernels::InitParams initParams = kernels::InitParams();
      initParams.latDat = mLatDat;
      initParams.lbmParams = &mParams;

      initParams.siteCount = mLatDat->GetInnerCollisionCount(0) + mLatDat->GetInterCollisionCount(0);
      mMidFluidCollision->Reset(&initParams);

      initParams.siteCount = mLatDat->GetInnerCollisionCount(1) + mLatDat->GetInterCollisionCount(1);
      mWallCollision->Reset(&initParams);

      initParams.siteCount = mLatDat->GetInnerCollisionCount(2) + mLatDat->GetInterCollisionCount(2);
      initParams.boundaryObject = mInletValues;
      mInletCollision->Reset(&initParams);

      initParams.siteCount = mLatDat->GetInnerCollisionCount(3) + mLatDat->GetInterCollisionCount(3);
      initParams.boundaryObject = mOutletValues;
      mOutletCollision->Reset(&initParams);

      initParams.siteCount = mLatDat->GetInnerCollisionCount(4) + mLatDat->GetInterCollisionCount(4);
      initParams.boundaryObject = mInletValues;
      mInletWallCollision->Reset(&initParams);

      initParams.siteCount = mLatDat->GetInnerCollisionCount(5) + mLatDat->GetInterCollisionCount(5);
      initParams.boundaryObject = mOutletValues;
      mOutletWallCollision->Reset(&initParams);
    }

    template<class LatticeType>
    LBM<LatticeType>::~LBM()
    {
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

    template<class LatticeType>
    void LBM<LatticeType>::ReadParameters()
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

      mParams.StressType = mSimConfig->StressType;
    }

    template<class LatticeType>
    void LBM<LatticeType>::WriteConfigParallel(hemelb::lb::Stability const stability,
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

      MPI_Datatype viewType = MpiDataType<char>();
      MPI_File_set_view(lOutputFile, 0, viewType, viewType, &lReadMode[0], MPI_INFO_NULL);

      topology::NetworkTopology* netTop = topology::NetworkTopology::Instance();

      if (netTop->IsCurrentProcTheIOProc())
      {
        // Write the header according to format detailed in snapshot.h
        char lBuffer[io::formats::snapshot::HeaderLength];
        io::writers::xdr::XdrMemWriter lWriter = io::writers::xdr::XdrMemWriter(lBuffer,
                                                                                io::formats::snapshot::HeaderLength);
        lWriter << (unsigned int) io::formats::HemeLbMagicNumber << (unsigned int) io::formats::snapshot::MagicNumber
            << (unsigned int) io::formats::snapshot::VersionNumber;
        lWriter << (unsigned int) io::formats::snapshot::HeaderLength;
        lWriter << stability;
        lWriter << mLatDat->GetVoxelSize();
        lWriter << mLatDat->GetOrigin().x << mLatDat->GetOrigin().x << mLatDat->GetOrigin().z;
        lWriter << (int) mLatDat->GetGlobalSiteMins().x << (int) mLatDat->GetGlobalSiteMins().y
            << (int) mLatDat->GetGlobalSiteMins().z;
        lWriter << (int) mLatDat->GetGlobalSiteMaxes().x << (int) mLatDat->GetGlobalSiteMaxes().y
            << (int) mLatDat->GetGlobalSiteMaxes().z;
        lWriter << (int) mLatDat->GetTotalFluidSites();

        MPI_File_write(lOutputFile, lBuffer, io::formats::snapshot::HeaderLength, MpiDataType(lBuffer[0]), &lStatus);
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
        lLocalSitesInitialOffset += io::formats::snapshot::VoxelRecordLength * mLatDat->GetFluidSiteCountOnProc(ii);
      }

      MPI_File_set_view(lOutputFile, lLocalSitesInitialOffset, viewType, viewType, &lReadMode[0], MPI_INFO_NULL);

      site_t lLocalWriteLength = io::formats::snapshot::VoxelRecordLength
          * mLatDat->GetFluidSiteCountOnProc(netTop->GetLocalRank());
      char * lFluidSiteBuffer = new char[lLocalWriteLength];
      hemelb::io::writers::xdr::XdrMemWriter lWriter =
          hemelb::io::writers::xdr::XdrMemWriter(lFluidSiteBuffer, (unsigned int) lLocalWriteLength);

      /* The following loops scan over every single macrocell (block). If
       * the block is non-empty, it scans the sites within that block. If the
       * site is fluid and present on the current task, it calculates the
       * flow field and encodes it to the local buffer.
       */
      for (geometry::BlockTraverser blockTrav(*mLatDat); blockTrav.CurrentLocationValid(); blockTrav.TraverseOne())
      {
        const geometry::Block& block = blockTrav.GetCurrentBlockData();

        if (block.IsEmpty())
        {
          continue;
        }

        for (geometry::SiteTraverser siteTrav = blockTrav.GetSiteTraverser(); siteTrav.CurrentLocationValid();
            siteTrav.TraverseOne())
        {
          if (netTop->GetLocalRank() != block.GetProcessorRankForSite(siteTrav.GetCurrentIndex()))
          {
            continue;
          }

          site_t my_site_id = block.GetLocalContiguousIndexForSite(siteTrav.GetCurrentIndex());

          /* Skip over solid sites. */
          if (my_site_id & BIG_NUMBER3)
            continue;

          distribn_t density, vx, vy, vz, f_eq[LatticeType::NUMVECTORS], f_neq[LatticeType::NUMVECTORS], stress,
              pressure;

          geometry::Site site = mLatDat->GetSite(my_site_id);

          if (site.GetSiteType() == geometry::FLUID_TYPE && !site.IsEdge())
          {
            LatticeType::CalculateDensityVelocityFEq(site.GetFOld(), density, vx, vy, vz, f_eq);

            for (unsigned int l = 0; l < LatticeType::NUMVECTORS; l++)
            {
              f_neq[l] = site.GetFOld()[l] - f_eq[l];
            }

          }
          else
          { // not FLUID_TYPE
            CalculateBC(site.GetFOld(), site.GetSiteType(), site.GetBoundaryId(), &density, &vx, &vy, &vz, f_neq);
          }

          if (mParams.StressType == hemelb::lb::ShearStress)
          {
            if (!site.IsEdge())
            {
              stress = -1.0;
            }
            else
            {
              LatticeType::CalculateShearStress(density,
                                                f_neq,
                                                site.GetWallNormal(),
                                                stress,
                                                mParams.GetStressParameter());
            }
          }
          else
          {
            LatticeType::CalculateVonMisesStress(f_neq, stress, mParams.GetStressParameter());
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

          const util::Vector3D<site_t> relativeSiteCoords = mLatDat->GetGlobalCoords(blockTrav.GetCurrentLocation(),
                                                                                     siteTrav.GetCurrentLocation())
              - siteMins;

          lWriter << (int) relativeSiteCoords.x << (int) relativeSiteCoords.y << (int) relativeSiteCoords.z;

          lWriter << float(pressure) << float(vx) << float(vy) << float(vz) << float(stress);
        }
      }

      if (netTop->GetProcessorCount() == 1)
      {
        // On hector, romio doesn't like to write_all to a single-machine communicator.
        // So we do a simple write.
        MPI_File_write(lOutputFile,
                       lFluidSiteBuffer,
                       (int) lLocalWriteLength,
                       MpiDataType(lFluidSiteBuffer[0]),
                       &lStatus);
      }
      else
      {
        // Hand the buffers over to MPIO to write to the file.
        MPI_File_write_all(lOutputFile,
                           lFluidSiteBuffer,
                           (int) lLocalWriteLength,
                           MpiDataType(lFluidSiteBuffer[0]),
                           &lStatus);
      }
      MPI_File_close(&lOutputFile);

      delete[] lFluidSiteBuffer;
    }

    // Calculate the BCs for each boundary site type and the
    // non-equilibrium distribution functions.
    template<class LatticeType>
    void LBM<LatticeType>::CalculateBC(distribn_t f[],
                                       hemelb::geometry::SiteType const iSiteType,
                                       unsigned int const iBoundaryId,
                                       distribn_t *density,
                                       distribn_t *vx,
                                       distribn_t *vy,
                                       distribn_t *vz,
                                       distribn_t f_neq[]) const
    {
      distribn_t dummy_density;

      for (unsigned int l = 0; l < LatticeType::NUMVECTORS; l++)
      {
        f_neq[l] = f[l];
      }

      // If you look at where this function is called from, having siteType == fluid type actually
      // means it also  must be an edge (pure fluid case is caught before this function is called).
      // UGH.
      if (iSiteType == hemelb::geometry::FLUID_TYPE)
      {
        LatticeType::CalculateDensityAndVelocity(f, *density, *vx, *vy, *vz);
      }
      else
      {
        if (iSiteType == hemelb::geometry::INLET_TYPE)
        {
          *density = mInletValues->GetBoundaryDensity(iBoundaryId);
        }
        else
        {
          *density = mOutletValues->GetBoundaryDensity(iBoundaryId);
        }

        LatticeType::CalculateDensityAndVelocity(f, dummy_density, *vx, *vy, *vz);
        LatticeType::CalculateFeq(*density, *vx, *vy, *vz, f);

      }
      for (unsigned int l = 0; l < LatticeType::NUMVECTORS; l++)
      {
        f_neq[l] -= f[l];
      }

    }

    template<class LatticeType>
    void LBM<LatticeType>::ReadVisParameters()
    {
      distribn_t density_min = std::numeric_limits < distribn_t > ::max();
      distribn_t density_max = std::numeric_limits < distribn_t > ::min();

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

#endif /* HEMELB_LB_LB_HPP */
