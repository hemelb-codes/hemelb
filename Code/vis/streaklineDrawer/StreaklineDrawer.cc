#include <math.h>
#include <vector>
#include <iostream>

#include "geometry/BlockTraverser.h"
#include "geometry/SiteTraverser.h"
#include "util/utilityFunctions.h"
#include "vis/streaklineDrawer/StreaklineDrawer.h"
#include "vis/Control.h"
#include "vis/streaklineDrawer/StreakPixel.h"
#include "vis/XYCoordinates.h"

namespace hemelb
{
  namespace vis
  {
    namespace streaklinedrawer
    {
      // TODO the streaker needs to be fixed. The lines drawn differ depending on how many proccesors
      // are used to do the visualisation. This is a bug.

      // Constructor, populating fields from lattice data objects.
      StreaklineDrawer::StreaklineDrawer(const geometry::LatticeData& iLatDat,
                                         const Screen& iScreen,
                                         const Viewpoint& iViewpoint, //
                                         const VisSettings& iVisSettings) :
        latDat(iLatDat), particleManager(neighbouringProcessors), mScreen(iScreen),
            velocityField(neighbouringProcessors), mViewpoint(iViewpoint),
            mVisSettings(iVisSettings)
      {

        site_t n = 0;

        topology::NetworkTopology* netTop = topology::NetworkTopology::Instance();

        velocityField.BuildVelocityField(iLatDat, this);

        procs = netTop->GetProcessorCount();
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
      void StreaklineDrawer::StreakLines(unsigned long time_steps,
                                         unsigned long time_steps_per_cycle)
      {
        // Set the particle creation period to be every time step, unless there are >=10000
        // timesteps per cycle
        unsigned int particle_creation_period =
            util::NumericalFunctions::max<unsigned int>(1, (unsigned int) (time_steps_per_cycle
                / 5000));

        int timestepsBetweenStreaklinesRounded = (int) (0.5F + (float) time_steps_per_cycle
            / mVisSettings.streaklines_per_pulsatile_period);

        if ((float) (time_steps % timestepsBetweenStreaklinesRounded)
            <= (mVisSettings.streakline_length / 100.0F) * ((float) time_steps_per_cycle
                / mVisSettings.streaklines_per_pulsatile_period) && time_steps
            % particle_creation_period == 0)
        {
          createSeedParticles();
        }

        bool debugThisIt = time_steps == 500 || time_steps == 1000;

        velocityField.InvalidateAllCalculatedVelocities();

        UpdateVelocityFieldForCommunicatedSites();

        CommunicateSiteIds();

        CommunicateVelocities();

        UpdateVelocityFieldForAllParticles();

        particleManager.ProcessParticleMovement();

        particleManager.CommunicateParticles(latDat, procs, velocityField);

        if (debugThisIt)
        {
          particleManager.PrintParticleCount();
        }
      }

      // Render the streaklines
      PixelSet<StreakPixel>* StreaklineDrawer::Render()
      {
        int pixels_x = mScreen.GetPixelsX();
        int pixels_y = mScreen.GetPixelsY();

        const std::vector<Particle>& particles = particleManager.GetParticles();

        PixelSet<StreakPixel>* set = GetUnusedPixelSet();
        set->Clear();

        for (unsigned int n = 0; n < particles.size(); n++)
        {
          util::Vector3D<float> p1;
          p1.x = particles[n].x - float (latDat.GetXSiteCount() >> 1);
          p1.y = particles[n].y - float (latDat.GetYSiteCount() >> 1);
          p1.z = particles[n].z - float (latDat.GetZSiteCount() >> 1);

          util::Vector3D<float> p2 = mViewpoint.Project(p1);

          XYCoordinates<int> x = mScreen.TransformScreenToPixelCoordinates<int> (XYCoordinates<
              float> (p2.x, p2.y));

          if (! (x.x < 0 || x.x >= pixels_x || x.y < 0 || x.y >= pixels_y))
          {
            StreakPixel pixel(x.x, x.y, particles[n].vel, p2.z, particles[n].inletID);
            set->AddPixel(pixel);
          }
        }

        return set;
      }

      // Create seed particles to begin the streaklines.
      void StreaklineDrawer::createSeedParticles()
      {
        for (unsigned int n = 0; n < particleSeeds.size(); n++)
        {
          particleManager.AddParticle(particleSeeds[n]);
        }
      }

      void StreaklineDrawer::UpdateVelocityFieldForAllParticles()
      {
        std::vector<Particle>& particles = particleManager.GetParticles();

        // Note that we iterate through the array back to front because at the end of this
        // loop we delete a particle. Iterating back to front ensures that we end up
        // visiting each particle.
        for (int n = (int) particles.size() - 1; n >= 0; --n)
        {
          float v[2][2][2][3];
          int is_interior;

          velocityField.GetVelocityFieldAroundPoint(util::Vector3D<site_t>((site_t) particles[n].x,
                                                                           (site_t) particles[n].y,
                                                                           (site_t) particles[n].z),
                                                    latDat,
                                                    v);

          float interp_v[3];

          velocityField.InterpolateVelocityForPoint(particles[n].x,
                                                    particles[n].y,
                                                    particles[n].z,
                                                    v,
                                                    interp_v);

          float vel = interp_v[0] * interp_v[0] + interp_v[1] * interp_v[1] + interp_v[2]
              * interp_v[2];

          if (vel > 1.0F)
          {
            particles[n].vel = 1.0F;
            particles[n].vx = interp_v[0] / sqrtf(vel);
            particles[n].vy = interp_v[1] / sqrtf(vel);
            particles[n].vz = interp_v[2] / sqrtf(vel);

          }
          else if (vel > 1.0e-8)
          {
            particles[n].vel = sqrtf(vel);
            particles[n].vx = interp_v[0];
            particles[n].vy = interp_v[1];
            particles[n].vz = interp_v[2];

          }
          else
          {
            particleManager.DeleteParticle(n);
          }
        }
      }

      void StreaklineDrawer::UpdateVelocityFieldForCommunicatedSites()
      {
        for (std::map<proc_t, NeighbouringProcessor>::iterator neighProc =
            neighbouringProcessors.begin(); neighProc != neighbouringProcessors.end(); neighProc++)
        {
          NeighbouringProcessor& proc = (*neighProc).second;

          for (site_t sendingVelocityIndex = 0; sendingVelocityIndex
              < proc.GetSitesToReceiveCount(); sendingVelocityIndex++)
          {
            const util::Vector3D<site_t>& siteCoords =
                proc.GetSiteCoordsBeingReceived(sendingVelocityIndex);

            float dummy1[2][2][2][3];

            velocityField.GetVelocityFieldAroundPoint(siteCoords, latDat, dummy1);
          }
        }
      }

      // Communicate site ids to other processors.
      void StreaklineDrawer::CommunicateSiteIds()
      {
        for (std::map<proc_t, NeighbouringProcessor>::iterator proc =
            neighbouringProcessors.begin(); proc != neighbouringProcessors.end(); proc++)
        {
          (*proc).second.ExchangeSiteIds();
        }
        for (std::map<proc_t, NeighbouringProcessor>::iterator proc =
            neighbouringProcessors.begin(); proc != neighbouringProcessors.end(); proc++)
        {
          (*proc).second.WaitForSiteIdExchange();
        }
      }

      // Communicate velocities to other processors.
      void StreaklineDrawer::CommunicateVelocities()
      {
        for (std::map<proc_t, NeighbouringProcessor>::iterator proc =
            neighbouringProcessors.begin(); proc != neighbouringProcessors.end(); proc++)
        {
          NeighbouringProcessor& neighbourProc = (*proc).second;

          for (site_t n = 0; n < neighbourProc.GetSitesToReceiveCount(); n++)
          {
            const util::Vector3D<site_t> siteCoords = neighbourProc.GetSiteCoordsBeingReceived(n);

            const VelocitySiteData* velocityDataForSite =
                velocityField.GetVelocitySiteData(latDat, siteCoords);

            neighbourProc.SetVelocityFieldToSend(n, velocityDataForSite != NULL
              ? util::Vector3D<float>(velocityDataForSite->vx,
                                      velocityDataForSite->vy,
                                      velocityDataForSite->vz)
              : util::Vector3D<float>(0., 0., 0.));
          }
        }

        for (std::map<proc_t, NeighbouringProcessor>::iterator proc =
            neighbouringProcessors.begin(); proc != neighbouringProcessors.end(); proc++)
        {
          (*proc).second.ExchangeVelocitiesForRequestedSites();
        }
        for (std::map<proc_t, NeighbouringProcessor>::iterator proc =
            neighbouringProcessors.begin(); proc != neighbouringProcessors.end(); proc++)
        {
          (*proc).second.WaitForVelocityExchange();
        }

        for (std::map<proc_t, NeighbouringProcessor>::iterator proc =
            neighbouringProcessors.begin(); proc != neighbouringProcessors.end(); proc++)
        {
          NeighbouringProcessor& neighbourProc = (*proc).second;

          if (neighbourProc.GetSendingSiteCount() == 0)
            continue;

          for (size_t n = 0; n < neighbourProc.GetSendingSiteCount(); n++)
          {
            util::Vector3D<site_t>& coords = neighbourProc.GetSiteCoorinates(n);

            VelocitySiteData* vel_site_data_p = velocityField.GetVelocitySiteData(latDat, coords);

            if (vel_site_data_p != NULL)
            {
              const util::Vector3D<float>& velocity = neighbourProc.GetReceivedVelocity(n);
              vel_site_data_p->vx = velocity.x;
              vel_site_data_p->vy = velocity.y;
              vel_site_data_p->vz = velocity.z;
            }
          }
        }

        for (std::map<proc_t, NeighbouringProcessor>::iterator proc =
            neighbouringProcessors.begin(); proc != neighbouringProcessors.end(); proc++)
        {
          (*proc).second.ClearSendingSites();
        }
      }

    }
  }
}
