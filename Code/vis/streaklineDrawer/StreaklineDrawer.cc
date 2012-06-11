#include <cmath>
#include <vector>
#include <iostream>

#include "geometry/BlockTraverser.h"
#include "geometry/SiteTraverser.h"
#include "util/utilityFunctions.h"
#include "vis/Control.h"
#include "vis/streaklineDrawer/StreaklineDrawer.h"
#include "vis/streaklineDrawer/StreakPixel.h"
#include "vis/XYCoordinates.h"

namespace hemelb
{
  namespace vis
  {
    namespace streaklinedrawer
    {
      // Constructor, populating fields from lattice data objects.
      StreaklineDrawer::StreaklineDrawer(const geometry::LatticeData& iLatDat,
                                         const Screen& iScreen,
                                         const Viewpoint& iViewpoint,
                                         const VisSettings& iVisSettings,
                                         const lb::MacroscopicPropertyCache& propertyCache) :
          latDat(iLatDat), screen(iScreen), viewpoint(iViewpoint), visSettings(iVisSettings), propertyCache(propertyCache), particleManager(neighbouringProcessors), velocityField(neighbouringProcessors,
                                                                                                                                                                                   propertyCache)
      {
        velocityField.BuildVelocityField(iLatDat);
        ChooseSeedParticles();
      }

      // Destructor
      StreaklineDrawer::~StreaklineDrawer()
      {
      }

      // Reset the streakline drawer.
      void StreaklineDrawer::Restart()
      {
        particleManager.DeleteAll();
      }

      // Create streakline particles and move them.
      void StreaklineDrawer::ProgressStreaklines(unsigned long time_steps, unsigned long total_time_steps)
      {
        // Set the particle creation period to be every time step, unless there are >=10000
        // timesteps
        unsigned int particle_creation_period =
            util::NumericalFunctions::max<unsigned int>(1, (unsigned int) (total_time_steps / 5000));

        int timestepsBetweenStreaklinesRounded = (int) (0.5F
            + (float) total_time_steps / visSettings.streaklines_per_simulation);

        if ((float) (time_steps % timestepsBetweenStreaklinesRounded)
            <= (visSettings.streakline_length / 100.0F)
                * ((float) total_time_steps / visSettings.streaklines_per_simulation)
            && time_steps % particle_creation_period == 0)
        {
          CreateParticlesFromSeeds();
        }

        velocityField.InvalidateAllCalculatedVelocities();

        // Decide which sites we will need velocity data for in order to move particles
        WorkOutVelocityDataNeededForParticles();

        // Communicate this with other processors.
        CommunicateSiteIds();

        // Recalculate velocities at the sites requested from this rank.
        UpdateVelocityFieldForCommunicatedSites();

        // Communicate the velocities back to the sites that requested them.
        CommunicateVelocities();

        // Update our local velocity field is sufficient to update all particles and prune
        // those that are stationary.
        UpdateVelocityFieldForAllParticlesAndPrune();

        // Process the particles' movement.
        particleManager.ProcessParticleMovement();

        // Communicate any particles that have crossed into the territory of another rank.
        particleManager.CommunicateParticles(latDat, velocityField);
      }

      // Render the streaklines
      PixelSet<StreakPixel>* StreaklineDrawer::Render()
      {
        int pixels_x = screen.GetPixelsX();
        int pixels_y = screen.GetPixelsY();

        const std::vector<Particle>& particles = particleManager.GetParticles();

        PixelSet<StreakPixel>* set = GetUnusedPixelSet();
        set->Clear();

        for (unsigned int n = 0; n < particles.size(); n++)
        {
          util::Vector3D<float> p1 = particles[n].position - util::Vector3D<float>(latDat.GetSiteDimensions() / 2);

          util::Vector3D<float> p2 = viewpoint.Project(p1);

          XYCoordinates<int> x = screen.TransformScreenToPixelCoordinates<int>(XYCoordinates<float>(p2.x, p2.y));

          if (! (x.x < 0 || x.x >= pixels_x || x.y < 0 || x.y >= pixels_y))
          {
            StreakPixel pixel(x.x, x.y, particles[n].vel, p2.z, particles[n].inletID);
            set->AddPixel(pixel);
          }
        }

        return set;
      }

      void StreaklineDrawer::ChooseSeedParticles()
      {
        site_t inlet_sites = 0;

        geometry::BlockTraverser blockTraverser(latDat);
        do
        {
          const geometry::Block& block = blockTraverser.GetCurrentBlockData();

          if (block.IsEmpty())
          {
            continue;
          }

          geometry::SiteTraverser siteTraverser(latDat);
          do
          {
            if (topology::NetworkTopology::Instance()->GetLocalRank()
                != block.GetProcessorRankForSite(siteTraverser.GetCurrentIndex()))
            {
              continue;
            }

            const geometry::ConstSite site =
                latDat.GetSite(block.GetLocalContiguousIndexForSite(siteTraverser.GetCurrentIndex()));

            // if the lattice site is not an inlet
            if (site.GetSiteType() != geometry::INLET_TYPE)
            {
              continue;
            }
            ++inlet_sites;

            // TODO this is a problem on multiple cores.
            if (inlet_sites % 50 != 0)
            {
              continue;
            }

            particleSeeds.push_back(Particle(static_cast<float>(blockTraverser.GetX() * blockTraverser.GetBlockSize()
                                                 + siteTraverser.GetX()),
                                             static_cast<float>(blockTraverser.GetY() * blockTraverser.GetBlockSize()
                                                 + siteTraverser.GetY()),
                                             static_cast<float>(blockTraverser.GetZ() * blockTraverser.GetBlockSize()
                                                 + siteTraverser.GetZ()),
                                             site.GetBoundaryId()));

          }
          while (siteTraverser.TraverseOne());
        }
        while (blockTraverser.TraverseOne());
      }

      // Create seed particles to begin the streaklines.
      void StreaklineDrawer::CreateParticlesFromSeeds()
      {
        for (unsigned int n = 0; n < particleSeeds.size(); n++)
        {
          particleManager.AddParticle(particleSeeds[n]);
        }
      }

      void StreaklineDrawer::WorkOutVelocityDataNeededForParticles()
      {
        std::vector<Particle> &particles = particleManager.GetParticles();

        // Note that we iterate through the array back to front because at the end of this
        // loop we delete a particle. Iterating back to front ensures that we end up
        // visiting each particle.
        for (int n = (int) particles.size() - 1; n >= 0; --n)
        {
          Particle& particle = particles[n];

          for (int unitGridI = 0; unitGridI <= 1; ++unitGridI)
          {
            site_t neighbourI = (site_t) particle.position.x + unitGridI;

            for (int unitGridJ = 0; unitGridJ <= 1; ++unitGridJ)
            {
              site_t neighbourJ = (site_t) particle.position.y + unitGridJ;

              for (int unitGridK = 0; unitGridK <= 1; ++unitGridK)
              {
                site_t neighbourK = (site_t) particle.position.z + unitGridK;

                proc_t sourceProcessor;

                if (velocityField.NeededFromNeighbour(util::Vector3D<site_t>(neighbourI, neighbourJ, neighbourK),
                                                      latDat,
                                                      &sourceProcessor))
                {
                  neighbouringProcessors[sourceProcessor].AddSiteToRequestVelocityDataFor(neighbourI,
                                                                                          neighbourJ,
                                                                                          neighbourK);
                }
              }
            }
          }
        }
      }

      void StreaklineDrawer::UpdateVelocityFieldForAllParticlesAndPrune()
      {
        std::vector<Particle> &particles = particleManager.GetParticles();

        // Note that we iterate through the array back to front because at the end of this
        // loop we delete a particle. Iterating back to front ensures that we end up
        // visiting each particle.
        for (int n = (int) particles.size() - 1; n >= 0; --n)
        {
          util::Vector3D<float> localVelocityField[2][2][2];

          velocityField.GetVelocityFieldAroundPoint(util::Vector3D<site_t>(particles[n].position),
                                                    latDat,
                                                    localVelocityField);

          util::Vector3D<float> interp_v = velocityField.InterpolateVelocityForPoint(particles[n].position,
                                                                                     localVelocityField);

          float vel = interp_v.Dot(interp_v);

          if (vel > 1.0F)
          {
            particles[n].vel = 1.0F;
            particles[n].velocity = interp_v * float(1.0 / sqrtf(vel));
          }
          else if (vel > 1.0e-8)
          {
            particles[n].vel = sqrtf(vel);
            particles[n].velocity = interp_v;
          }
          else
          {
            particleManager.DeleteParticle(n);
          }
        }
      }

      void StreaklineDrawer::UpdateVelocityFieldForCommunicatedSites()
      {
        for (std::map<proc_t, NeighbouringProcessor>::const_iterator neighProc = neighbouringProcessors.begin();
            neighProc != neighbouringProcessors.end(); ++neighProc)
        {
          const NeighbouringProcessor& proc = (*neighProc).second;

          for (site_t sendingVelocityIndex = 0; sendingVelocityIndex < proc.GetNumberOfSitesRequestedByNeighbour();
              sendingVelocityIndex++)
          {
            const util::Vector3D<site_t>& siteCoords =
                proc.GetSiteCoordsBeingRequestedByNeighbour(sendingVelocityIndex);

            velocityField.UpdateLocalField(siteCoords, latDat);
          }
        }
      }

      // Communicate site ids to other processors.
      void StreaklineDrawer::CommunicateSiteIds()
      {
        net::Net net;

        for (std::map<proc_t, NeighbouringProcessor>::iterator proc = neighbouringProcessors.begin();
            proc != neighbouringProcessors.end(); ++proc)
        {
          (*proc).second.ExchangeSiteIdCounts(net);
        }

        net.Receive();
        net.Send();
        net.Wait();

        for (std::map<proc_t, NeighbouringProcessor>::iterator proc = neighbouringProcessors.begin();
            proc != neighbouringProcessors.end(); proc++)
        {
          (*proc).second.ExchangeSiteIds(net);
        }

        net.Receive();
        net.Send();
        net.Wait();
      }

      // Communicate velocities to other processors.
      void StreaklineDrawer::CommunicateVelocities()
      {
        net::Net net;

        for (std::map<proc_t, NeighbouringProcessor>::iterator proc = neighbouringProcessors.begin();
            proc != neighbouringProcessors.end(); ++proc)
        {
          (*proc).second.ExchangeVelocitiesForRequestedSites(net);
        }

        for (std::map<proc_t, NeighbouringProcessor>::iterator proc = neighbouringProcessors.begin();
            proc != neighbouringProcessors.end(); ++proc)
        {
          NeighbouringProcessor& neighbourProc = (*proc).second;

          for (site_t n = 0; n < neighbourProc.GetNumberOfSitesRequestedByNeighbour(); ++n)
          {
            const util::Vector3D<site_t> siteCoords = neighbourProc.GetSiteCoordsBeingRequestedByNeighbour(n);

            const VelocitySiteData* velocityDataForSite = velocityField.GetVelocitySiteData(latDat, siteCoords);

            neighbourProc.SetVelocityFieldToSend(n, velocityDataForSite->velocity);
          }
        }

        net.Receive();
        net.Send();
        net.Wait();

        for (std::map<proc_t, NeighbouringProcessor>::const_iterator proc = neighbouringProcessors.begin();
            proc != neighbouringProcessors.end(); ++proc)
        {
          const NeighbouringProcessor& neighbourProc = (*proc).second;

          for (site_t n = 0; n < neighbourProc.GetNumberOfSitesRequestedByThisCore(); n++)
          {
            const util::Vector3D<site_t> &coords = neighbourProc.GetSendingSiteCoorinates(n);
            velocityField.GetVelocitySiteData(latDat, coords)->velocity = neighbourProc.GetReceivedVelocityField(n);
          }
        }

        for (std::map<proc_t, NeighbouringProcessor>::iterator proc = neighbouringProcessors.begin();
            proc != neighbouringProcessors.end(); ++proc)
        {
          (*proc).second.ClearListOfRequestedSites();
        }
      }

    }
  }
}
