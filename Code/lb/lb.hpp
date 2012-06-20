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
    lb::MacroscopicPropertyCache& LBM<LatticeType>::GetPropertyCache()
    {
      return propertyCache;
    }

    template<class LatticeType>
    LBM<LatticeType>::LBM(configuration::SimConfig *iSimulationConfig,
                          net::Net* net,
                          geometry::LatticeData* latDat,
                          SimulationState* simState,
                          reporting::Timers &atimings,
                          geometry::neighbouring::NeighbouringDataManager *neighbouringDataManager) :
        mSimConfig(iSimulationConfig), mNet(net), mLatDat(latDat), mState(simState), mParams(mState->GetTimeStepLength(),
                                                                                             latDat->GetVoxelSize()), timings(atimings), propertyCache(*simState,
                                                                                                                                                       *latDat), neighbouringDataManager(neighbouringDataManager)
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
      initParams.neighbouringDataManager=neighbouringDataManager;

      initParams.siteCount = mLatDat->GetMidDomainCollisionCount(0) + mLatDat->GetDomainEdgeCollisionCount(0);
      mMidFluidCollision = new tMidFluidCollision(initParams);

      initParams.siteCount = mLatDat->GetMidDomainCollisionCount(1) + mLatDat->GetDomainEdgeCollisionCount(1);
      mWallCollision = new tWallCollision(initParams);

      initParams.siteCount = mLatDat->GetMidDomainCollisionCount(2) + mLatDat->GetDomainEdgeCollisionCount(2);
      initParams.boundaryObject = mInletValues;
      mInletCollision = new tInletOutletCollision(initParams);

      initParams.siteCount = mLatDat->GetMidDomainCollisionCount(3) + mLatDat->GetDomainEdgeCollisionCount(3);
      initParams.boundaryObject = mOutletValues;
      mOutletCollision = new tInletOutletCollision(initParams);

      initParams.siteCount = mLatDat->GetMidDomainCollisionCount(4) + mLatDat->GetDomainEdgeCollisionCount(4);
      initParams.boundaryObject = mInletValues;
      mInletWallCollision = new tInletOutletWallCollision(initParams);

      initParams.siteCount = mLatDat->GetMidDomainCollisionCount(5) + mLatDat->GetDomainEdgeCollisionCount(5);
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

      for (int i = 0; i < OutletCount(); i++)
      {
        density += mOutletValues->GetDensityMin(i);
      }

      density /= OutletCount();

      for (site_t i = 0; i < mLatDat->GetLocalFluidSiteCount(); i++)
      {
        distribn_t f_eq[LatticeType::NUMVECTORS];

        LatticeType::CalculateFeq(density, 0.0, 0.0, 0.0, f_eq);

        geometry::Site site = mLatDat->GetSite(i);

        distribn_t* f_old_p = site.GetFOld<LatticeType>();
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
      timings[hemelb::reporting::Timers::lb].Start();

      // Delegate to the lattice data object to post the asynchronous sends and receives
      // (via the Net object).
      // NOTE that this doesn't actually *perform* the sends and receives, it asks the Net
      // to include them in the ISends and IRecvs that happen later.
      mLatDat->SendAndReceive(mNet);

      timings[hemelb::reporting::Timers::lb].Stop();
    }

    template<class LatticeType>
    void LBM<LatticeType>::PreSend()
    {
      timings[hemelb::reporting::Timers::lb].Start();
      timings[hemelb::reporting::Timers::lb_calc].Start();

      /**
       * In the PreSend phase, we do LB on all the sites that need to have results sent to
       * neighbouring ranks ('domainEdge' sites). In site id terms, this means we start at the
       * end of the sites whose neighbours all lie on this rank ('midDomain'), then progress
       * through the sites of each type in turn.
       */
      site_t offset = mLatDat->GetMidDomainSiteCount();

      StreamAndCollide(mMidFluidCollision, offset, mLatDat->GetDomainEdgeCollisionCount(0));
      offset += mLatDat->GetDomainEdgeCollisionCount(0);

      StreamAndCollide(mWallCollision, offset, mLatDat->GetDomainEdgeCollisionCount(1));
      offset += mLatDat->GetDomainEdgeCollisionCount(1);

      mInletValues->FinishReceive();
      StreamAndCollide(mInletCollision, offset, mLatDat->GetDomainEdgeCollisionCount(2));
      offset += mLatDat->GetDomainEdgeCollisionCount(2);

      mOutletValues->FinishReceive();
      StreamAndCollide(mOutletCollision, offset, mLatDat->GetDomainEdgeCollisionCount(3));
      offset += mLatDat->GetDomainEdgeCollisionCount(3);

      StreamAndCollide(mInletWallCollision, offset, mLatDat->GetDomainEdgeCollisionCount(4));
      offset += mLatDat->GetDomainEdgeCollisionCount(4);

      StreamAndCollide(mOutletWallCollision, offset, mLatDat->GetDomainEdgeCollisionCount(5));

      timings[hemelb::reporting::Timers::lb].Stop();
    }

    template<class LatticeType>
    void LBM<LatticeType>::PreReceive()
    {
      timings[hemelb::reporting::Timers::lb].Start();

      /**
       * In the PreReceive phase, we perform LB for all the sites whose neighbours lie on this
       * rank ('midDomain' rather than 'domainEdge' sites). Ideally this phase is the longest bit (maximising time for the asynchronous sends
       * and receives to complete).
       *
       * In site id terms, this means starting at the first site and progressing through the
       * midDomain sites, one type at a time.
       */
      site_t offset = 0;

      StreamAndCollide(mMidFluidCollision, offset, mLatDat->GetMidDomainCollisionCount(0));
      offset += mLatDat->GetMidDomainCollisionCount(0);

      StreamAndCollide(mWallCollision, offset, mLatDat->GetMidDomainCollisionCount(1));
      offset += mLatDat->GetMidDomainCollisionCount(1);

      StreamAndCollide(mInletCollision, offset, mLatDat->GetMidDomainCollisionCount(2));
      offset += mLatDat->GetMidDomainCollisionCount(2);

      StreamAndCollide(mOutletCollision, offset, mLatDat->GetMidDomainCollisionCount(3));
      offset += mLatDat->GetMidDomainCollisionCount(3);

      StreamAndCollide(mInletWallCollision, offset, mLatDat->GetMidDomainCollisionCount(4));
      offset += mLatDat->GetMidDomainCollisionCount(4);

      StreamAndCollide(mOutletWallCollision, offset, mLatDat->GetMidDomainCollisionCount(5));

      timings[hemelb::reporting::Timers::lb_calc].Stop();
      timings[hemelb::reporting::Timers::lb].Stop();
    }

    template<class LatticeType>
    void LBM<LatticeType>::PostReceive()
    {
      timings[hemelb::reporting::Timers::lb].Start();

      // Copy the distribution functions received from the neighbouring
      // processors into the destination buffer "f_new".
      // This is done here, after receiving the sent distributions from neighbours.
      mLatDat->CopyReceived();

      // Do any cleanup steps necessary on boundary nodes
      site_t offset = 0;

      timings[hemelb::reporting::Timers::lb_calc].Start();

      //TODO yup, this is horrible. If you read this, please improve the following code.
      PostStep(mMidFluidCollision, offset, mLatDat->GetDomainEdgeCollisionCount(0));
      offset += mLatDat->GetDomainEdgeCollisionCount(0);

      PostStep(mWallCollision, offset, mLatDat->GetDomainEdgeCollisionCount(1));
      offset += mLatDat->GetDomainEdgeCollisionCount(1);

      PostStep(mInletCollision, offset, mLatDat->GetDomainEdgeCollisionCount(2));
      offset += mLatDat->GetDomainEdgeCollisionCount(2);

      PostStep(mOutletCollision, offset, mLatDat->GetDomainEdgeCollisionCount(3));
      offset += mLatDat->GetDomainEdgeCollisionCount(3);

      PostStep(mInletWallCollision, offset, mLatDat->GetDomainEdgeCollisionCount(4));
      offset += mLatDat->GetDomainEdgeCollisionCount(4);

      PostStep(mOutletWallCollision, offset, mLatDat->GetDomainEdgeCollisionCount(5));
      offset += mLatDat->GetDomainEdgeCollisionCount(5);

      PostStep(mMidFluidCollision, offset, mLatDat->GetMidDomainCollisionCount(0));
      offset += mLatDat->GetMidDomainCollisionCount(0);

      PostStep(mWallCollision, offset, mLatDat->GetMidDomainCollisionCount(1));
      offset += mLatDat->GetMidDomainCollisionCount(1);

      PostStep(mInletCollision, offset, mLatDat->GetMidDomainCollisionCount(2));
      offset += mLatDat->GetMidDomainCollisionCount(2);

      PostStep(mOutletCollision, offset, mLatDat->GetMidDomainCollisionCount(3));
      offset += mLatDat->GetMidDomainCollisionCount(3);

      PostStep(mInletWallCollision, offset, mLatDat->GetMidDomainCollisionCount(4));
      offset += mLatDat->GetMidDomainCollisionCount(4);

      PostStep(mOutletWallCollision, offset, mLatDat->GetMidDomainCollisionCount(5));

      timings[hemelb::reporting::Timers::lb_calc].Stop();
      timings[hemelb::reporting::Timers::lb].Stop();
    }

    template<class LatticeType>
    void LBM<LatticeType>::EndIteration()
    {
      timings[hemelb::reporting::Timers::lb].Start();

      // Swap f_old and f_new ready for the next timestep.
      mLatDat->SwapOldAndNew();

      timings[hemelb::reporting::Timers::lb].Stop();
    }

    // In the case of instability, this function restart the simulation
    // with twice as many time steps per period and update the parameters
    // that depends on this change.
    template<class LatticeType>
    void LBM<LatticeType>::Reset()
    {
      mState->DoubleTimeResolution();

      mParams.Update(mState->GetTimeStepLength(), mLatDat->GetVoxelSize());

      SetInitialConditions();

      kernels::InitParams initParams = kernels::InitParams();
      initParams.latDat = mLatDat;
      initParams.lbmParams = &mParams;

      initParams.siteCount = mLatDat->GetMidDomainCollisionCount(0) + mLatDat->GetDomainEdgeCollisionCount(0);
      mMidFluidCollision->Reset(&initParams);

      initParams.siteCount = mLatDat->GetMidDomainCollisionCount(1) + mLatDat->GetDomainEdgeCollisionCount(1);
      mWallCollision->Reset(&initParams);

      initParams.siteCount = mLatDat->GetMidDomainCollisionCount(2) + mLatDat->GetDomainEdgeCollisionCount(2);
      initParams.boundaryObject = mInletValues;
      mInletCollision->Reset(&initParams);

      initParams.siteCount = mLatDat->GetMidDomainCollisionCount(3) + mLatDat->GetDomainEdgeCollisionCount(3);
      initParams.boundaryObject = mOutletValues;
      mOutletCollision->Reset(&initParams);

      initParams.siteCount = mLatDat->GetMidDomainCollisionCount(4) + mLatDat->GetDomainEdgeCollisionCount(4);
      initParams.boundaryObject = mInletValues;
      mInletWallCollision->Reset(&initParams);

      initParams.siteCount = mLatDat->GetMidDomainCollisionCount(5) + mLatDat->GetDomainEdgeCollisionCount(5);
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
    }

    template<class LatticeType>
    void LBM<LatticeType>::ReadParameters()
    {
      std::vector<lb::boundaries::iolets::InOutLet*> inlets = mSimConfig->GetInlets();
      std::vector<lb::boundaries::iolets::InOutLet*> outlets = mSimConfig->GetOutlets();
      inletCount = inlets.size();
      outletCount = outlets.size();
      mParams.StressType = mSimConfig->GetStressType();
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
        lWriter << (int) stability;
        lWriter << (double) mLatDat->GetVoxelSize();
        lWriter << (double) mLatDat->GetOrigin().x << (double) mLatDat->GetOrigin().y
            << (double) mLatDat->GetOrigin().z;
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

          /* Skip over solid sites. */
          if (block.SiteIsSolid(siteTrav.GetCurrentIndex()))
          {
            continue;
          }

          site_t my_site_id = block.GetLocalContiguousIndexForSite(siteTrav.GetCurrentIndex());

          distribn_t stress;

          geometry::Site site = mLatDat->GetSite(my_site_id);

          /// @todo #111 It should possible to compute shear stress in the whole domain, not only at the boundary
          if (mParams.StressType == hemelb::lb::ShearStress && !site.IsEdge())
          {
            /**
             *  @todo #111 This is a pretty meaningless way of saying that you cannot compute stress for this site.
             *  The -1 value will be later on translated to physical units and will end up being different values
             *  for different runs depending on grid and simulation parameters.
             */
            stress = -1.0;
          }
          else if (mParams.StressType)
          {
            /// @todo #138, wall normals don't seem to be initialised properly, they appear to be [nan,nan,nan]
            stress = propertyCache.shearStressCache.Get(my_site_id);
          }
          else
          {
            stress = propertyCache.vonMisesStressCache.Get(my_site_id);
          }

          // conversion from lattice to physical units
          distribn_t pressure = mUnits->ConvertPressureToPhysicalUnits(propertyCache.densityCache.Get(my_site_id)
              * Cs2);

          const util::Vector3D<distribn_t>& velocity = propertyCache.velocityCache.Get(my_site_id);

          distribn_t vx = mUnits->ConvertVelocityToPhysicalUnits(velocity.x);
          distribn_t vy = mUnits->ConvertVelocityToPhysicalUnits(velocity.y);
          distribn_t vz = mUnits->ConvertVelocityToPhysicalUnits(velocity.z);

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

    template<class LatticeType>
    void LBM<LatticeType>::ReadVisParameters()
    {
      distribn_t density_min = std::numeric_limits<distribn_t>::max();
      distribn_t density_max = std::numeric_limits<distribn_t>::min();

      distribn_t velocity_max = mUnits->ConvertVelocityToLatticeUnits(mSimConfig->GetMaximumVelocity());
      distribn_t stress_max = mUnits->ConvertStressToLatticeUnits(mSimConfig->GetMaximumStress());

      for (int i = 0; i < InletCount(); i++)
      {
        density_min = util::NumericalFunctions::min(density_min, mInletValues->GetDensityMin(i));
        density_max = util::NumericalFunctions::max(density_max, mInletValues->GetDensityMax(i));
      }
      for (int i = 0; i < OutletCount(); i++)
      {
        density_min = util::NumericalFunctions::min(density_min, mOutletValues->GetDensityMin(i));
        density_max = util::NumericalFunctions::max(density_max, mOutletValues->GetDensityMax(i));
      }

      distribn_t lDensity_threshold_min = density_min;
      distribn_t lDensity_threshold_minmax_inv = 1.0F / (density_max - density_min);
      distribn_t lVelocity_threshold_max_inv = 1.0F / velocity_max;
      distribn_t lStress_threshold_max_inv = 1.0F / stress_max;

      mVisControl->SetSomeParams(mSimConfig->GetVisualisationBrightness(),
                                 lDensity_threshold_min,
                                 lDensity_threshold_minmax_inv,
                                 lVelocity_threshold_max_inv,
                                 lStress_threshold_max_inv);
    }
  }
}

#endif /* HEMELB_LB_LB_HPP */
