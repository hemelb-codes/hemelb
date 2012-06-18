#include "colloids/Particle.h"
#include "geometry/LatticeData.h"
#include "log/Logger.h"

namespace hemelb
{
  namespace colloids
  {
    Particle::Particle(const geometry::LatticeData* const latDatLBM,
                       io::xml::XmlAbstractionLayer& xml,
                       lb::MacroscopicPropertyCache& propertyCache) :
      PersistedParticle(xml),
      latDatLBM(latDatLBM),
      propertyCache(propertyCache)
    {
    }

    const void Particle::UpdatePosition()
    {
      if (log::Logger::ShouldDisplay<log::Debug>())
        log::Logger::Log<log::Debug, log::OnePerCore>(
          "In colloids::Particle::UpdatePosition, id: %i,\nposition: {%g,%g,%g}\nvelocity: {%g,%g,%g}\nbodyForces: {%g,%g,%g}\n",
          topology::NetworkTopology::Instance()->GetLocalRank(),
          particleId, globalPosition.x, globalPosition.y, globalPosition.z,
          velocity.x, velocity.y, velocity.z, bodyForces.x, bodyForces.y, bodyForces.z);

      globalPosition += velocity + bodyForces;

      if (log::Logger::ShouldDisplay<log::Debug>())
        log::Logger::Log<log::Debug, log::OnePerCore>(
          "In colloids::Particle::UpdatePosition, id: %i, position is now: {%g,%g,%g}\n",
          topology::NetworkTopology::Instance()->GetLocalRank(),
          particleId, globalPosition.x, globalPosition.y, globalPosition.z);
    }

    const void Particle::CalculateBodyForces()
    {
      /** todo: CalculateBodyForces
       *    - create new BodyForces::GetBodyForcesAtGlobalPosition method
       *    - set bodyForces property to value from the BodyForces object
       *    - will require communicaiton to transmit remote contributions
       */
      bodyForces *= 0.0;
    }

    const void Particle::CalculateFeedbackForces() const
    {
      /** todo: CalculateFeedbackForces
       *    - calculate the feedback force on each neighbour lattice site
       *    - set feedback force values into the macroscopic cache object
       *    - will require communication to transmit remote contributions
       */
    }

    /** modified dirac delta function according to Peskin */
    double diracOperation(util::Vector3D<double> relativePosition)
    {
      double delta = 1.0;
      for (int xyz=0;xyz<3;xyz++)
      {
        double r = relativePosition[xyz];
        double rmod = fabs(r);

        if (rmod <= 1.0)
          delta *= 0.125*(3.0 - 2.0*rmod + sqrt(1.0 + 4.0*rmod - 4.0*rmod*rmod));
        else if (rmod <= 2.0)
          delta *= 0.125*(5.0 - 2.0*rmod - sqrt(-7.0 + 12.0*rmod  - 4.0*rmod*rmod));
        else
          delta = 0.0;

      }
      return delta;
    }

    const void Particle::InterpolateFluidVelocity()
    {
      /** InterpolateFluidVelocity
       *    For each local neighbour lattice site
       *    - get velocity for the site from the macroscopic cache object
       *    - calculate site's contribution to the velocity interpolation
       *    - increment particle velocity property with this contribution
       *    - will require communication to transmit remote contributions
       */

      if (log::Logger::ShouldDisplay<log::Debug>())
        log::Logger::Log<log::Debug, log::OnePerCore>(
          "In colloids::Particle::InterpolateFluidVelocity, id: %i, position: {%g,%g,%g}\n",
          particleId, globalPosition.x, globalPosition.y, globalPosition.z);

      velocity *= 0.0;
      // determine the global coordinates of the next neighbour site:
      // nested loop - x, y, z directions semi-open interval [-2, +2)
      for (int x = ((int)globalPosition.x)-1;
               x < ((int)globalPosition.x)+2; x++)
        for (int y = ((int)globalPosition.y)-1;
                 y < ((int)globalPosition.y)+2; y++)
          for (int z = ((int)globalPosition.z)-1;
                   z < ((int)globalPosition.z)+2; z++)
          {
            const util::Vector3D<site_t> siteGlobalPosition(x, y, z);

            // convert the global coordinates of the site into a local site index
            proc_t procId;
            site_t siteId;
            bool isValid = latDatLBM->GetContiguousSiteId(siteGlobalPosition, procId, siteId);
            bool isLocal = (procId == topology::NetworkTopology::Instance()->GetLocalRank());

            if (log::Logger::ShouldDisplay<log::Debug>())
            {
              util::Vector3D<site_t> blockCoords, localSiteCoords;
              latDatLBM->GetBlockAndLocalSiteCoords(siteGlobalPosition, blockCoords, localSiteCoords);

              log::Logger::Log<log::Debug, log::OnePerCore>(
                "In colloids::Particle::InterpolateFluidVelocity, particleId: %i,\nsiteGlobalPosition: {%i,%i,%i}\nblockCoords: {%i,%i,%i} - %s\nlocalSiteCoords: {%i,%i,%i} - %s\nprocId: %i, siteId: %i\n",
                particleId,
                siteGlobalPosition.x, siteGlobalPosition.y, siteGlobalPosition.z,
                blockCoords.x, blockCoords.y, blockCoords.z,
                !latDatLBM->IsValidBlock(blockCoords) ? "INVALID" :
                (latDatLBM->GetBlock(latDatLBM->GetBlockIdFromBlockCoords(blockCoords)).IsEmpty() ?
                "SOLID" : "valid"),
                localSiteCoords.x, localSiteCoords.y, localSiteCoords.z,
                !latDatLBM->IsValidLatticeSite(localSiteCoords) ? "INVALID" :
                (!isValid ? "SOLID" : (!isLocal ? "REMOTE" : "local fluid")),
                procId, siteId);
            }

            /** TODO: implement communication for remote fluid sites */
            /** TODO: implement boundary conditions for invalid/solid sites */
            if (!isValid || !isLocal)
              continue;

            // read value of velocity for site index from macroscopic cache
            util::Vector3D<double> siteFluidVelocity = propertyCache.velocityCache.Get(siteId);

            // calculate term of the interpolation sum
            util::Vector3D<double> relativePosition(siteGlobalPosition);
            relativePosition -= globalPosition;
            util::Vector3D<double> partialInterpolation = siteFluidVelocity *
              diracOperation(relativePosition);

            // accumulate each term of the interpolation
            velocity += partialInterpolation;

            if (log::Logger::ShouldDisplay<log::Debug>())
              log::Logger::Log<log::Debug, log::OnePerCore>(
                "In colloids::Particle::InterpolateFluidVelocity, particleId: %i, siteIndex: %i, fluidVelocity: {%g,%g,%g}, partialInterpolation: {%g,%g,%g}, velocitySoFar: {%g,%g,%g}\n",
                particleId, siteId, siteFluidVelocity.x, siteFluidVelocity.y, siteFluidVelocity.z,
                partialInterpolation.x, partialInterpolation.y, partialInterpolation.z,
                velocity.x, velocity.y, velocity.z);
          }
      propertyCache.velocityCache.SetRefreshFlag();
    }

  }
}
