#include <vector>
#include <math.h>
#include <limits>

#include "log/Logger.h"
#include "util/utilityFunctions.h"
#include "vis/Control.h"
#include "vis/rayTracer/RayTracer.h"
#include "vis/GlyphDrawer.h"

#include "io/XdrFileWriter.h"

namespace hemelb
		
{
  namespace vis
  {
    Control::Control(lb::StressTypes iStressType,
                     net::Net* net,
                     lb::SimulationState* simState,
                     geometry::LatticeData* iLatDat) :
      net::PhasedBroadcastIrregular<true, 2, 0, false, true>(net, simState, SPREADFACTOR),
          mLatDat(iLatDat)
    {
      timeSpent = 0.0;

      mVisSettings.mStressType = iStressType;

      this->vis = new Vis;

      //sites_x etc are globals declared in net.h
      vis->half_dim[0] = 0.5F * float (iLatDat->GetXSiteCount());
      vis->half_dim[1] = 0.5F * float (iLatDat->GetYSiteCount());
      vis->half_dim[2] = 0.5F * float (iLatDat->GetZSiteCount());

      vis->system_size = 2.F * fmaxf(vis->half_dim[0], fmaxf(vis->half_dim[1], vis->half_dim[2]));

      mVisSettings.mouse_x = -1;
      mVisSettings.mouse_y = -1;

      initLayers();
    }

    void Control::initLayers()
    {
      // We don't have all the minima / maxima on one core, so we have to gather them.
      // NOTE this only happens once, during initialisation, otherwise it would be
      // totally unforgivable.
      site_t block_min_x = std::numeric_limits<site_t>::max();
      site_t block_min_y = std::numeric_limits<site_t>::max();
      site_t block_min_z = std::numeric_limits<site_t>::max();
      site_t block_max_x = std::numeric_limits<site_t>::min();
      site_t block_max_y = std::numeric_limits<site_t>::min();
      site_t block_max_z = std::numeric_limits<site_t>::min();

      site_t n = -1;

      for (site_t i = 0; i < mLatDat->GetXBlockCount(); i++)
      {
        for (site_t j = 0; j < mLatDat->GetYBlockCount(); j++)
        {
          for (site_t k = 0; k < mLatDat->GetZBlockCount(); k++)
          {
            n++;

            geometry::LatticeData::BlockData * lBlock = mLatDat->GetBlock(n);
            if (lBlock->ProcessorRankForEachBlockSite == NULL)
            {
              continue;
            }

            block_min_x = util::NumericalFunctions::min(block_min_x, i);
            block_min_y = util::NumericalFunctions::min(block_min_y, j);
            block_min_z = util::NumericalFunctions::min(block_min_z, k);
            block_max_x = util::NumericalFunctions::max(block_max_x, i);
            block_max_y = util::NumericalFunctions::max(block_max_y, j);
            block_max_z = util::NumericalFunctions::max(block_max_z, k);
          }
        }
      }

      site_t mins[3], maxes[3];
      site_t localMins[3], localMaxes[3];

      localMins[0] = block_min_x;
      localMins[1] = block_min_y;
      localMins[2] = block_min_z;

      localMaxes[0] = block_max_x;
      localMaxes[1] = block_max_y;
      localMaxes[2] = block_max_z;

      MPI_Allreduce(localMins, mins, 3, MpiDataType<site_t> (), MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(localMaxes, maxes, 3, MpiDataType<site_t> (), MPI_MAX, MPI_COMM_WORLD);

      mVisSettings.ctr_x = 0.5F * (float) (mLatDat->GetBlockSize() * (mins[0] + maxes[0]));
      mVisSettings.ctr_y = 0.5F * (float) (mLatDat->GetBlockSize() * (mins[1] + maxes[1]));
      mVisSettings.ctr_z = 0.5F * (float) (mLatDat->GetBlockSize() * (mins[2] + maxes[2]));

      myRayTracer = new raytracer::RayTracer(mLatDat, &mDomainStats, &mScreen, &mViewpoint, &mVisSettings);
      myGlypher = new GlyphDrawer(mLatDat, &mScreen, &mDomainStats, &mViewpoint, &mVisSettings);

#ifndef NO_STREAKLINES
      myStreaker = new StreaklineDrawer(mLatDat, &mScreen, &mViewpoint, &mVisSettings);
#endif
      // Note that rtInit does stuff to this->ctr_x (because this has
      // to be global)
      mVisSettings.ctr_x -= vis->half_dim[0];
      mVisSettings.ctr_y -= vis->half_dim[1];
      mVisSettings.ctr_z -= vis->half_dim[2];
    }

    void Control::SetProjection(const int &iPixels_x,
                                const int &iPixels_y,
                                const float &iLocal_ctr_x,
                                const float &iLocal_ctr_y,
                                const float &iLocal_ctr_z,
                                const float &iLongitude,
                                const float &iLatitude,
                                const float &iZoom)
    {
      float rad = 5.F * vis->system_size;
      float dist = 0.5F * rad;

      Vector3D<float> centre = Vector3D<float>( iLocal_ctr_x, iLocal_ctr_y, iLocal_ctr_z );

      mViewpoint.SetViewpointPosition(iLongitude * (float) DEG_TO_RAD, iLatitude
          * (float) DEG_TO_RAD, centre, rad, dist);

      mScreen.Set( (0.5F * vis->system_size) / iZoom,
                   (0.5F * vis->system_size) / iZoom,
                  iPixels_x,
                  iPixels_y,
                  rad,
                  &mViewpoint);
    }

    void Control::RegisterSite(site_t i, distribn_t density, distribn_t velocity, distribn_t stress)
    {
      myRayTracer->UpdateClusterVoxel(i, density, velocity, stress);
    }

    void Control::SetSomeParams(const float iBrightness,
                                const distribn_t iDensityThresholdMin,
                                const distribn_t iDensityThresholdMinMaxInv,
                                const distribn_t iVelocityThresholdMaxInv,
                                const distribn_t iStressThresholdMaxInv)
    {
      mVisSettings.brightness = iBrightness;
      mDomainStats.density_threshold_min = iDensityThresholdMin;

      mDomainStats.density_threshold_minmax_inv = iDensityThresholdMinMaxInv;
      mDomainStats.velocity_threshold_max_inv = iVelocityThresholdMaxInv;
      mDomainStats.stress_threshold_max_inv = iStressThresholdMaxInv;
    }

    void Control::UpdateImageSize(int pixels_x, int pixels_y)
    {
      mScreen.Resize(pixels_x, pixels_y);
    }

    void Control::Render()
    {
      mScreen.Reset();

      myRayTracer->Render();

      if (mVisSettings.mode == VisSettings::ISOSURFACESANDGLYPHS)
      {
        myGlypher->Render();
      }
#ifndef NO_STREAKLINES
      if (mVisSettings.mStressType == lb::ShearStress || mVisSettings.mode
          == VisSettings::WALLANDSTREAKLINES)
      {
        myStreaker->render(mLatDat);
      }
#endif

      log::Logger::Log<log::Debug, log::OnePerCore>("Rendering.");
    }

    void Control::InitialAction(unsigned long startIteration)
    {
      timeSpent -= util::myClock();

      Render();

      log::Logger::Log<log::Debug, log::OnePerCore>("Render stored for phased imaging.");

      ScreenPixels* pix;

      // If we don't have any in the buffer, create a new ScreenPixels object.
      if (pixelsBuffer.empty())
      {
        pix = mScreen.SwapBuffers(new ScreenPixels());
      }
      // Otherwise use a ScreenPixels object from the buffer.
      else
      {
        ScreenPixels* buff = pixelsBuffer.top();
        pixelsBuffer.pop();
        pix = mScreen.SwapBuffers(buff);
      }

      resultsByStartIt.insert(std::pair<unsigned long, ScreenPixels*>(startIteration, pix));

      timeSpent += util::myClock();
    }

    void Control::ProgressFromChildren(unsigned long startIteration, unsigned long splayNumber)
    {
      timeSpent -= util::myClock();

      if (splayNumber == 0)
      {
        unsigned int* childNumbers[SPREADFACTOR];
        unsigned int counts[SPREADFACTOR];
        for (unsigned int ii = 0; ii < SPREADFACTOR; ++ii)
        {
          recvBuffers[ii].Reset();
          childNumbers[ii] = recvBuffers[ii].GetStoredPixelCountPtr();
          counts[ii] = 1;
        }

        log::Logger::Log<log::Debug, log::OnePerCore>("Receiving child image pixel count.");

        ReceiveFromChildren<unsigned int> (childNumbers, counts);
      }
      else if (splayNumber == 1)
      {
        ColPixel* childData[SPREADFACTOR];
        unsigned int counts[SPREADFACTOR];
        for (unsigned int ii = 0; ii < SPREADFACTOR; ++ii)
        {
          childData[ii] = recvBuffers[ii].GetPixelArray();
          counts[ii] = recvBuffers[ii].GetStoredPixelCount();
        }

        log::Logger::Log<log::Debug, log::OnePerCore>("Receiving child image pixel data.");

        ReceiveFromChildren<ColPixel> (childData, counts);
      }

      timeSpent += util::myClock();
    }

    void Control::ProgressToParent(unsigned long startIteration, unsigned long splayNumber)
    {
      timeSpent -= util::myClock();

      ScreenPixels* pixels = resultsByStartIt[startIteration];
      if (splayNumber == 0)
      {
        log::Logger::Log<log::Debug, log::OnePerCore>("Sending pixel count.");

        SendToParent<unsigned int> (pixels->GetStoredPixelCountPtr(), 1);
      }
      else if (splayNumber == 1)
      {
        log::Logger::Log<log::Debug, log::OnePerCore>("Sending pixel data.");

        SendToParent<ColPixel> (pixels->GetPixelArray(), pixels->GetStoredPixelCount());
      }

      timeSpent += util::myClock();
    }

    void Control::PostReceiveFromChildren(unsigned long startIteration, unsigned long splayNumber)
    {
      timeSpent -= util::myClock();

      if (splayNumber == 1)
      {
        ScreenPixels* pixels = resultsByStartIt[startIteration];

        log::Logger::Log<log::Debug, log::OnePerCore>("Combining in child pixel data.");

        for (unsigned int child = 0; child < SPREADFACTOR; ++child)
        {
          pixels->FoldIn(&recvBuffers[child], &mVisSettings);
        }
      }

      timeSpent += util::myClock();
    }

    bool Control::IsRendering() const
    {
      return IsInitialAction() || IsInstantBroadcast();
    }

    void Control::ClearOut(unsigned long startIt)
    {
      timeSpent -= util::myClock();

      bool found;

      do
      {
        found = false;

        if (resultsByStartIt.size() > 0)
        {
          mapType::iterator it = resultsByStartIt.begin();
          if (it->first <= startIt)
          {
            log::Logger::Log<log::Debug, log::OnePerCore>("Clearing out image cache from it %lu",
                                                          it->first);

            it->second->Reset();
            pixelsBuffer.push(it->second);
            resultsByStartIt.erase(it);
            found = true;
          }
        }
      }
      while (found);

      timeSpent += util::myClock();
    }

    const ScreenPixels* Control::GetResult(unsigned long startIt)
    {
      log::Logger::Log<log::Debug, log::OnePerCore>("Getting image results from it %lu", startIt);

      if (resultsByStartIt.count(startIt) != 0)
      {
        return resultsByStartIt[startIt];
      }
      else
      {
        return NULL;
      }
    }

    void Control::InstantBroadcast(unsigned long startIteration)
    {
      timeSpent -= util::myClock();

      log::Logger::Log<log::Debug, log::OnePerCore>("Performing instant imaging.");

      Render();

      // Status object for MPI comms.
      MPI_Status status;

      /*
       * We do several iterations.
       *
       * On the first, every even proc passes data to the odd proc below, where it is merged.
       * On the second, the difference is two, so proc 3 passes to 1, 7 to 5, 11 to 9 etc.
       * On the third the differenec is four, so proc 5 passes to 1, 13 to 9 etc.
       * .
       * .
       * .
       *
       * This continues until all data is passed back to processor one, which passes it to proc 0.
       */
      topology::NetworkTopology* netTop = topology::NetworkTopology::Instance();

      // Start with a difference in rank of 1, doubling every time.
      for (proc_t deltaRank = 1; deltaRank < netTop->GetProcessorCount(); deltaRank <<= 1)
      {
        // The receiving proc is all the ranks that are 1 modulo (deltaRank * 2)
        for (proc_t receivingProc = 1; receivingProc < (netTop->GetProcessorCount() - deltaRank); receivingProc
            += deltaRank << 1)
        {
          proc_t sendingProc = receivingProc + deltaRank;

          // If we're the sending proc, do the send.
          if (netTop->GetLocalRank() == sendingProc)
          {
            MPI_Send(mScreen.pixels->GetStoredPixelCountPtr(),
                     1,
                     MpiDataType(mScreen.pixels->GetStoredPixelCount()),
                     receivingProc,
                     20,
                     MPI_COMM_WORLD);

            if (mScreen.pixels->GetStoredPixelCount() > 0)
            {
              MPI_Send(mScreen.pixels->GetPixelArray(),
                       mScreen.pixels->GetStoredPixelCount(),
                       MpiDataType<ColPixel> (),
                       receivingProc,
                       20,
                       MPI_COMM_WORLD);
            }
          }

          // If we're the receiving proc, receive.
          else if (netTop->GetLocalRank() == receivingProc)
          {
            MPI_Recv(recvBuffers[0].GetStoredPixelCountPtr(),
                     1,
                     MpiDataType(recvBuffers[0].GetStoredPixelCount()),
                     sendingProc,
                     20,
                     MPI_COMM_WORLD,
                     &status);

            if (recvBuffers[0].GetStoredPixelCount() > 0)
            {
              MPI_Recv(recvBuffers[0].GetPixelArray(),
                       recvBuffers[0].GetStoredPixelCount(),
                       MpiDataType<ColPixel> (),
                       sendingProc,
                       20,
                       MPI_COMM_WORLD,
                       &status);

              mScreen.pixels->FoldIn(&recvBuffers[0], &mVisSettings);
            }
          }
        }
      }

      // Send the final image from proc 1 to 0.
      if (netTop->GetLocalRank() == 1)
      {
        MPI_Send(mScreen.pixels->GetStoredPixelCountPtr(),
                 1,
                 MpiDataType(mScreen.pixels->GetStoredPixelCount()),
                 0,
                 20,
                 MPI_COMM_WORLD);

        if (mScreen.pixels->GetStoredPixelCount() > 0)
        {
          MPI_Send(mScreen.pixels->GetPixelArray(),
                   mScreen.pixels->GetStoredPixelCount(),
                   MpiDataType<ColPixel> (),
                   0,
                   20,
                   MPI_COMM_WORLD);
        }

      }
      // Receive the final image on proc 0.
      else if (netTop->GetLocalRank() == 0)
      {
        MPI_Recv(recvBuffers[0].GetStoredPixelCountPtr(),
                 1,
                 MpiDataType(recvBuffers[0].GetStoredPixelCount()),
                 1,
                 20,
                 MPI_COMM_WORLD,
                 &status);

        if (recvBuffers[0].GetStoredPixelCount() > 0)
        {
          MPI_Recv(recvBuffers[0].GetPixelArray(),
                   recvBuffers[0].GetStoredPixelCount(),
                   MpiDataType<ColPixel> (),
                   1,
                   20,
                   MPI_COMM_WORLD,
                   &status);

          mScreen.pixels->FoldIn(&recvBuffers[0], &mVisSettings);
        }

        ScreenPixels* pix;

        // Create a new pixels object if we don't have any spare ones in the buffer.
        if (pixelsBuffer.empty())
        {
          pix = mScreen.SwapBuffers(new ScreenPixels());
        }
        // Use a pixels object from the buffer when there is one.
        else
        {
          ScreenPixels* newBuff = pixelsBuffer.top();
          pixelsBuffer.pop();
          pix = mScreen.SwapBuffers(newBuff);
        }

        resultsByStartIt.insert(std::pair<unsigned long, ScreenPixels*>(base::mSimState->GetTimeStepsPassed(),
                                                                        pix));

        log::Logger::Log<log::Debug, log::OnePerCore>("Inserting image at it %lu.",
                                                      base::mSimState->GetTimeStepsPassed());
      }

      timeSpent += util::myClock();
    }

    void Control::SetMouseParams(double iPhysicalPressure, double iPhysicalStress)
    {
      mVisSettings.mouse_pressure = iPhysicalPressure;
      mVisSettings.mouse_stress = iPhysicalStress;
    }

    bool Control::MouseIsOverPixel(float* density, float* stress)
    {
      if (mVisSettings.mouse_x < 0 || mVisSettings.mouse_y < 0)
      {
        return false;
      }

      return mScreen.MouseIsOverPixel(mVisSettings.mouse_x, mVisSettings.mouse_y, density, stress);
    }

    void Control::ProgressStreaklines(unsigned long time_step, unsigned long period)
    {
      timeSpent -= util::myClock();

      myStreaker ->StreakLines(time_step, period, mLatDat);

      timeSpent += util::myClock();
    }

    double Control::GetTimeSpent() const
    {
      return timeSpent;
    }

    void Control::Reset()
    {
      timeSpent = 0.0;

      log::Logger::Log<log::Debug, log::OnePerCore>("Resetting image controller.");

      myStreaker->Restart();

      base::Reset();

      ClearOut(base::mSimState->GetTimeStepsPassed() + 1);
    }

    Control::~Control()
    {
#ifndef NO_STREAKLINES
      delete myStreaker;
#endif

      delete vis;
      delete myGlypher;
      delete myRayTracer;

      // Clear out the ScreenPixels used still in the results buffer.
      for (std::map<unsigned long, ScreenPixels*>::iterator it = resultsByStartIt.begin(); it
          != resultsByStartIt.end(); it++)
      {
        delete it->second;
      }

      // Clear out the ScreenPixels objects from the pixelsBuffer.
      while (!pixelsBuffer.empty())
      {
        delete pixelsBuffer.top();
        pixelsBuffer.pop();
      }
    }

  } // namespace vis
}
// namespace hemelb
