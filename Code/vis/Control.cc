// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#include <vector>
#include <cmath>
#include <limits>

#include "log/Logger.h"
#include "util/utilityFunctions.h"
#include "vis/Control.h"
#include "vis/rayTracer/RayTracer.h"
#include "vis/GlyphDrawer.h"

#include "io/writers/xdr/XdrFileWriter.h"

namespace hemelb

{
  namespace vis
  {
    Control::Control(lb::StressTypes iStressType,
                     net::Net* netIn,
                     lb::SimulationState* simState,
                     const lb::MacroscopicPropertyCache& propertyCache,
                     geometry::LatticeData* iLatDat,
                     reporting::Timer &atimer) :
        net::PhasedBroadcastIrregular<true, 2, 0, false, true>(netIn, simState, SPREADFACTOR),
        net(netIn), propertyCache(propertyCache), latticeData(iLatDat), timer(atimer)
    {

      visSettings.mStressType = iStressType;

      this->vis = new Vis;

      //sites_x etc are globals declared in net.h
      vis->half_dim = util::Vector3D<float>(iLatDat->GetSiteDimensions()) * 0.5F;

      vis->system_size = 2.F * fmaxf(vis->half_dim[0], fmaxf(vis->half_dim[1], vis->half_dim[2]));

      visSettings.mouse_x = -1;
      visSettings.mouse_y = -1;

      initLayers();
    }

    void Control::initLayers()
    {
      // We don't have all the minima / maxima on one core, so we have to gather them.
      // NOTE this only happens once, during initialisation, otherwise it would be
      // totally unforgivable.
      site_t block_min_x = std::numeric_limits < site_t > ::max();
      site_t block_min_y = std::numeric_limits < site_t > ::max();
      site_t block_min_z = std::numeric_limits < site_t > ::max();
      site_t block_max_x = std::numeric_limits < site_t > ::min();
      site_t block_max_y = std::numeric_limits < site_t > ::min();
      site_t block_max_z = std::numeric_limits < site_t > ::min();

      for (geometry::BlockTraverser blockIt(*latticeData); blockIt.CurrentLocationValid(); blockIt.TraverseOne())
      {
        if (blockIt.GetCurrentBlockData().IsEmpty())
        {
          continue;
        }

        block_min_x = util::NumericalFunctions::min(block_min_x, blockIt.GetCurrentLocation().x);
        block_min_y = util::NumericalFunctions::min(block_min_y, blockIt.GetCurrentLocation().y);
        block_min_z = util::NumericalFunctions::min(block_min_z, blockIt.GetCurrentLocation().z);
        block_max_x = util::NumericalFunctions::max(block_max_x, blockIt.GetCurrentLocation().x);
        block_max_y = util::NumericalFunctions::max(block_max_y, blockIt.GetCurrentLocation().y);
        block_max_z = util::NumericalFunctions::max(block_max_z, blockIt.GetCurrentLocation().z);
      }

      site_t mins[3], maxes[3];
      site_t localMins[3], localMaxes[3];

      localMins[0] = block_min_x;
      localMins[1] = block_min_y;
      localMins[2] = block_min_z;

      localMaxes[0] = block_max_x;
      localMaxes[1] = block_max_y;
      localMaxes[2] = block_max_z;

      MPI_Allreduce(localMins, mins, 3, MpiDataType<site_t>(), MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(localMaxes, maxes, 3, MpiDataType<site_t>(), MPI_MAX, MPI_COMM_WORLD);

      visSettings.ctr_x = 0.5F * (float) (latticeData->GetBlockSize() * (mins[0] + maxes[0]));
      visSettings.ctr_y = 0.5F * (float) (latticeData->GetBlockSize() * (mins[1] + maxes[1]));
      visSettings.ctr_z = 0.5F * (float) (latticeData->GetBlockSize() * (mins[2] + maxes[2]));

      normalRayTracer =
          new raytracer::RayTracer<raytracer::ClusterWithWallNormals, raytracer::RayDataNormal>(latticeData,
                                                                                                &domainStats,
                                                                                                &screen,
                                                                                                &viewpoint,
                                                                                                &visSettings);

      myGlypher = new GlyphDrawer(latticeData, &screen, &domainStats, &viewpoint, &visSettings);

#ifndef NO_STREAKLINES
      myStreaker = new streaklinedrawer::StreaklineDrawer(*latticeData, screen, viewpoint, visSettings, propertyCache);
#else
      myStreaker = NULL;
#endif
      // Note that rtInit does stuff to this->ctr_x (because this has
      // to be global)
      visSettings.ctr_x -= vis->half_dim[0];
      visSettings.ctr_y -= vis->half_dim[1];
      visSettings.ctr_z -= vis->half_dim[2];
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

      //For now set the maximum draw distance to twice the radius;
      visSettings.maximumDrawDistance = 2.0F * rad;

      util::Vector3D<float> centre = util::Vector3D<float>(iLocal_ctr_x, iLocal_ctr_y, iLocal_ctr_z);

      viewpoint.SetViewpointPosition(iLongitude * (float) DEG_TO_RAD,
                                     iLatitude * (float) DEG_TO_RAD,
                                     centre,
                                     rad,
                                     dist);

      screen.Set( (0.5F * vis->system_size) / iZoom,
                 (0.5F * vis->system_size) / iZoom,
                 iPixels_x,
                 iPixels_y,
                 rad,
                 &viewpoint);
    }

    void Control::SetSomeParams(const float iBrightness,
                                const distribn_t iDensityThresholdMin,
                                const distribn_t iDensityThresholdMinMaxInv,
                                const distribn_t iVelocityThresholdMaxInv,
                                const distribn_t iStressThresholdMaxInv)
    {
      visSettings.brightness = iBrightness;
      domainStats.density_threshold_min = iDensityThresholdMin;

      domainStats.density_threshold_minmax_inv = iDensityThresholdMinMaxInv;
      domainStats.velocity_threshold_max_inv = iVelocityThresholdMaxInv;
      domainStats.stress_threshold_max_inv = iStressThresholdMaxInv;
    }

    void Control::UpdateImageSize(int pixels_x, int pixels_y)
    {
      screen.Resize(pixels_x, pixels_y);
    }

    void Control::Render(unsigned long startIteration)
    {
      log::Logger::Log<log::Debug, log::OnePerCore>("Rendering.");

      PixelSet<raytracer::RayDataNormal>* ray = normalRayTracer->Render(propertyCache);

      PixelSet<BasicPixel> *glyph = NULL;

      if (visSettings.mode == VisSettings::ISOSURFACESANDGLYPHS)
      {
        glyph = myGlypher->Render(propertyCache);
      }
      else
      {
        glyph = myGlypher->GetUnusedPixelSet();
        glyph->Clear();
      }

      PixelSet<streaklinedrawer::StreakPixel> *streak = NULL;

      if (myStreaker != NULL
          && (visSettings.mStressType == lb::ShearStress || visSettings.mode == VisSettings::WALLANDSTREAKLINES))
      {
        streak = myStreaker->Render();
      }

      localResultsByStartIt.insert(std::pair<unsigned long, Rendering>(startIteration, Rendering(glyph, ray, streak)));
    }

    void Control::InitialAction(unsigned long startIteration)
    {
      timer.Start();

      Render(startIteration);

      log::Logger::Log<log::Debug, log::OnePerCore>("Render stored for phased imaging.");

      timer.Stop();
    }

    void Control::WriteImage(io::writers::Writer* writer,
                             const PixelSet<ResultPixel>& imagePixels,
                             const DomainStats& domainStats,
                             const VisSettings& visSettings) const
    {
      *writer << (int) visSettings.mode;

      *writer << domainStats.physical_pressure_threshold_min << domainStats.physical_pressure_threshold_max
          << domainStats.physical_velocity_threshold_max << domainStats.physical_stress_threshold_max;

      *writer << screen.GetPixelsX();
      *writer << screen.GetPixelsY();
      *writer << (int) imagePixels.GetPixelCount();

      WritePixels(writer, imagePixels, domainStats, visSettings);
    }

    int Control::GetPixelsX() const
    {
      return screen.GetPixelsX();
    }

    int Control::GetPixelsY() const
    {
      return screen.GetPixelsY();
    }

    void Control::WritePixels(io::writers::Writer* writer,
                              const PixelSet<ResultPixel>& imagePixels,
                              const DomainStats& domainStats,
                              const VisSettings& visSettings) const
    {
      const int bits_per_char = sizeof(char) * 8;

      for (unsigned int i = 0; i < imagePixels.GetPixelCount(); i++)
      {
        const ResultPixel& pixel = imagePixels.GetPixels()[i];

        // Use a ray-tracer function to get the necessary pixel data.
        unsigned index;
        unsigned char rgb_data[12];

        pixel.WritePixel(&index, rgb_data, domainStats, visSettings);

        *writer << (uint32_t) index;

        unsigned pix_data[3];
        pix_data[0] = (rgb_data[0] << (3 * bits_per_char)) + (rgb_data[1] << (2 * bits_per_char))
            + (rgb_data[2] << bits_per_char) + rgb_data[3];

        pix_data[1] = (rgb_data[4] << (3 * bits_per_char)) + (rgb_data[5] << (2 * bits_per_char))
            + (rgb_data[6] << bits_per_char) + rgb_data[7];

        pix_data[2] = (rgb_data[8] << (3 * bits_per_char)) + (rgb_data[9] << (2 * bits_per_char))
            + (rgb_data[10] << bits_per_char) + rgb_data[11];

        for (int i = 0; i < 3; i++)
        {
          *writer << (uint32_t) pix_data[i];
        }
        *writer << io::writers::Writer::eol;
      }
    }

    void Control::ProgressFromChildren(unsigned long startIteration, unsigned long splayNumber)
    {
      timer.Start();

      if (splayNumber == 0)
      {
        for (unsigned int ii = 0; ii < GetChildren().size(); ++ii)
        {
          Rendering lRendering(myGlypher->GetUnusedPixelSet(),
                               normalRayTracer->GetUnusedPixelSet(),
                               myStreaker != NULL ?
                                 myStreaker->GetUnusedPixelSet() :
                                 NULL);

          lRendering.ReceivePixelCounts(net, GetChildren()[ii]);

          childrenResultsByStartIt.insert(std::pair<unsigned long, Rendering>(startIteration, lRendering));
        }

        log::Logger::Log<log::Debug, log::OnePerCore>("Receiving child image pixel count.");
      }
      else if (splayNumber == 1)
      {
        std::multimap<unsigned long, Rendering>::iterator renderings =
            childrenResultsByStartIt.equal_range(startIteration).first;

        for (unsigned int ii = 0; ii < GetChildren().size(); ++ii)
        {
          Rendering& received = (*renderings).second;

          log::Logger::Log<log::Debug, log::OnePerCore>("Receiving child image pixel data (from it %li).",
                                                        startIteration);

          received.ReceivePixelData(net, GetChildren()[ii]);

          renderings++;
        }
      }

      timer.Stop();
    }

    void Control::ProgressToParent(unsigned long startIteration, unsigned long splayNumber)
    {
      timer.Start();

      Rendering& rendering = (*localResultsByStartIt.find(startIteration)).second;
      if (splayNumber == 0)
      {
        log::Logger::Log<log::Debug, log::OnePerCore>("Sending pixel count (from it %li).", startIteration);

        rendering.SendPixelCounts(net, GetParent());
      }
      else if (splayNumber == 1)
      {
        log::Logger::Log<log::Debug, log::OnePerCore>("Sending pixel data (from it %li).", startIteration);

        rendering.SendPixelData(net, GetParent());
      }

      timer.Stop();
    }

    void Control::PostReceiveFromChildren(unsigned long startIteration, unsigned long splayNumber)
    {
      timer.Start();

      // The first time round, ensure that we have enough memory to receive the image data that
      // will come next time.
      if (splayNumber == 0)
      {

      }
      if (splayNumber == 1)
      {
        std::pair < std::multimap<unsigned long, Rendering>::iterator , std::multimap<unsigned long, Rendering>::iterator
            > its = childrenResultsByStartIt.equal_range(startIteration);

        Rendering local = (*localResultsByStartIt.find(startIteration)).second;

        std::multimap<unsigned long, Rendering>::iterator rendering = its.first;
        while (rendering != its.second)
        {
          local.Combine( (*rendering).second);

          (*rendering).second.ReleaseAll();

          rendering++;
        }

        if (its.first != its.second)
        {
          its.first++;
          childrenResultsByStartIt.erase(its.first, its.second);
        }

        log::Logger::Log<log::Debug, log::OnePerCore>("Combining in child pixel data.");
      }

      timer.Stop();
    }

    void Control::PostSendToParent(unsigned long startIteration, unsigned long splayNumber)
    {
      if (splayNumber == 1)
      {
        Rendering& rendering = (*localResultsByStartIt.find(startIteration)).second;
        rendering.ReleaseAll();

        localResultsByStartIt.erase(startIteration);
      }
    }

    bool Control::IsRendering() const
    {
      return IsInitialAction() || IsInstantBroadcast();
    }

    void Control::ClearOut(unsigned long startIt)
    {
      timer.Start();

      bool found;

      do
      {
        found = false;

        if (localResultsByStartIt.size() > 0)
        {
          mapType::iterator it = localResultsByStartIt.begin();
          if (it->first <= startIt)
          {
            log::Logger::Log<log::Debug, log::OnePerCore>("Clearing out image cache from it %lu", it->first);

            (*it).second.ReleaseAll();

            localResultsByStartIt.erase(it);
            found = true;
          }
        }
      }
      while (found);

      do
      {
        found = false;

        if (childrenResultsByStartIt.size() > 0)
        {
          multimapType::iterator it = childrenResultsByStartIt.begin();
          if ( (*it).first <= startIt)
          {
            log::Logger::Log<log::Debug, log::OnePerCore>("Clearing out image cache from it %lu", (*it).first);

            (*it).second.ReleaseAll();

            childrenResultsByStartIt.erase(it);
            found = true;
          }
        }
      }
      while (found);

      do
      {
        found = false;

        if (renderingsByStartIt.size() > 0)
        {
          std::multimap<unsigned long, PixelSet<ResultPixel>*>::iterator it = renderingsByStartIt.begin();
          if ( (*it).first <= startIt)
          {
            log::Logger::Log<log::Debug, log::OnePerCore>("Clearing out image cache from it %lu", (*it).first);

            (*it).second->Release();
            renderingsByStartIt.erase(it);
            found = true;
          }
        }
      }
      while (found);

      timer.Stop();
    }

    const PixelSet<ResultPixel>* Control::GetResult(unsigned long startIt)
    {
      log::Logger::Log<log::Debug, log::OnePerCore>("Getting image results from it %lu", startIt);

      if (renderingsByStartIt.count(startIt) != 0)
      {
        return (*renderingsByStartIt.find(startIt)).second;
      }

      if (localResultsByStartIt.count(startIt) != 0)
      {
        Rendering finalRender = (*localResultsByStartIt.find(startIt)).second;
        PixelSet<ResultPixel> *result = GetUnusedPixelSet();

        finalRender.PopulateResultSet(result);

        renderingsByStartIt.insert(std::pair<unsigned long, PixelSet<ResultPixel>*>(startIt, result));
        return result;
      }
      else
      {
        return NULL;
      }
    }

    void Control::InstantBroadcast(unsigned long startIteration)
    {
      timer.Start();

      log::Logger::Log<log::Debug, log::OnePerCore>("Performing instant imaging.");

      Render(startIteration);

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
      net::Net tempNet;

      Rendering* localBuffer = localResultsByStartIt.count(startIteration) > 0 ?
        & (*localResultsByStartIt.find(startIteration)).second :
        NULL;
      Rendering receiveBuffer(myGlypher->GetUnusedPixelSet(), normalRayTracer->GetUnusedPixelSet(), myStreaker == NULL ?
        NULL :
        myStreaker->GetUnusedPixelSet());

      // Start with a difference in rank of 1, doubling every time.
      for (proc_t deltaRank = 1; deltaRank < netTop->GetProcessorCount(); deltaRank <<= 1)
      {
        // The receiving proc is all the ranks that are 1 modulo (deltaRank * 2)
        for (proc_t receivingProc = 1; receivingProc < (netTop->GetProcessorCount() - deltaRank);
            receivingProc += deltaRank << 1)
        {
          proc_t sendingProc = receivingProc + deltaRank;

          // If we're the sending proc, do the send.
          if (netTop->GetLocalRank() == sendingProc)
          {
            localBuffer->SendPixelCounts(&tempNet, receivingProc);

            tempNet.Dispatch();

            localBuffer->SendPixelData(&tempNet, receivingProc);

            tempNet.Dispatch();
          }

          // If we're the receiving proc, receive.

          else if (netTop->GetLocalRank() == receivingProc)
          {
            receiveBuffer.ReceivePixelCounts(&tempNet, sendingProc);

            tempNet.Dispatch();

            receiveBuffer.ReceivePixelData(&tempNet, sendingProc);

            tempNet.Dispatch();

            localBuffer->Combine(receiveBuffer);
          }
        }
      }

      // Send the final image from proc 1 to 0.
      if (netTop->GetLocalRank() == 1)
      {
        localBuffer->SendPixelCounts(&tempNet, 0);

        tempNet.Dispatch();

        localBuffer->SendPixelData(&tempNet, 0);

        tempNet.Dispatch();
      }
      // Receive the final image on proc 0.

      else if (netTop->GetLocalRank() == 0)
      {
        receiveBuffer.ReceivePixelCounts(&tempNet, 1);

        tempNet.Dispatch();

        receiveBuffer.ReceivePixelData(&tempNet, 1);

        tempNet.Dispatch();

        localResultsByStartIt.erase(startIteration);
        localResultsByStartIt.insert(std::pair<unsigned long, Rendering>(startIteration, Rendering(receiveBuffer)));

        log::Logger::Log<log::Debug, log::OnePerCore>("Inserting image at it %lu.", startIteration);
      }

      if (netTop->GetLocalRank() != 0)
      {
        receiveBuffer.ReleaseAll();
      }

      timer.Stop();
    }

    void Control::SetMouseParams(double iPhysicalPressure, double iPhysicalStress)
    {
      visSettings.mouse_pressure = iPhysicalPressure;
      visSettings.mouse_stress = iPhysicalStress;
    }

    bool Control::MouseIsOverPixel(const PixelSet<ResultPixel>* result, float* density, float* stress)
    {
      if (visSettings.mouse_x < 0 || visSettings.mouse_y < 0)
      {
        return false;
      }

      const std::vector<ResultPixel>& screenPix = result->GetPixels();

      for (std::vector<ResultPixel>::const_iterator it = screenPix.begin(); it != screenPix.end(); ++it)
      {
        if ( (*it).GetRayPixel() != NULL && (*it).GetI() == visSettings.mouse_x && (*it).GetJ() == visSettings.mouse_y)
        {
          *density = (*it).GetRayPixel()->GetNearestDensity();
          *stress = (*it).GetRayPixel()->GetNearestStress();

          return true;
        }
      }

      return false;
    }

    void Control::ProgressStreaklines(unsigned long time_step, unsigned long period)
    {
      if (myStreaker != NULL)
      {
        timer.Start();
        myStreaker->ProgressStreaklines(time_step, period);
        timer.Stop();
      }
    }

    Control::~Control()
    {
      delete myStreaker;
      delete vis;
      delete myGlypher;
      delete normalRayTracer;
    }

  } // namespace vis
}
// namespace hemelb
