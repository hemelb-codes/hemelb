#include "colloids/Particle.h"
#include "geometry/LatticeData.h"
#include "log/Logger.h"
#include "units.h"

namespace hemelb
{
  namespace colloids
  {
    Particle::Particle(const geometry::LatticeData* const latDatLBM,
                       io::xml::XmlAbstractionLayer& xml,
                       lb::MacroscopicPropertyCache& propertyCache) :
      PersistedParticle(xml), latDatLBM(latDatLBM), propertyCache(propertyCache)
    {
      /** TODO: this conversion should be done in the xml abstraction layer */
      smallRadius_a0 /= latDatLBM->GetVoxelSize();
      largeRadius_ah /= latDatLBM->GetVoxelSize();
      globalPosition -= latDatLBM->GetOrigin();
      globalPosition /= latDatLBM->GetVoxelSize();

      log::Logger::Log<log::Trace, log::OnePerCore>("In colloids::Particle::ctor, id: %i, a0: %g, ah: %g, position: {%g,%g,%g}\n",
                                                    particleId,
                                                    smallRadius_a0,
                                                    largeRadius_ah,
                                                    globalPosition.x,
                                                    globalPosition.y,
                                                    globalPosition.z);
    }

    const void Particle::UpdatePosition()
    {
      log::Logger::Log<log::Trace, log::OnePerCore>("In colloids::Particle::UpdatePosition, id: %i,\nposition: {%g,%g,%g}\nvelocity: {%g,%g,%g}\nbodyForces: {%g,%g,%g}\n",
                                                    particleId,
                                                    globalPosition.x,
                                                    globalPosition.y,
                                                    globalPosition.z,
                                                    velocity.x,
                                                    velocity.y,
                                                    velocity.z,
                                                    bodyForces.x,
                                                    bodyForces.y,
                                                    bodyForces.z);

      globalPosition += velocity + bodyForces;

      log::Logger::Log<log::Trace, log::OnePerCore>("In colloids::Particle::UpdatePosition, id: %i, position is now: {%g,%g,%g}\n",
                                                    particleId,
                                                    globalPosition.x,
                                                    globalPosition.y,
                                                    globalPosition.z);
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
    //double diracOperation(util::Vector3D<double> relativePosition)
    DimensionlessQuantity diracOperation(LatticePosition relativePosition)
    {
      DimensionlessQuantity delta = 1.0;
      for (int xyz = 0; xyz < 3; xyz++)
      {
        LatticeDistance rmod = fabs(relativePosition[xyz]);

        if (rmod <= 1.0)
          delta *= 0.125 * (3.0 - 2.0 * rmod + sqrt(1.0 + 4.0 * rmod - 4.0 * rmod * rmod));
        else if (rmod <= 2.0)
          delta *= 0.125 * (5.0 - 2.0 * rmod - sqrt(-7.0 + 12.0 * rmod - 4.0 * rmod * rmod));
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

      log::Logger::Log<log::Trace, log::OnePerCore>("In colloids::Particle::InterpolateFluidVelocity, id: %i, position: {%g,%g,%g}\n",
                                                    particleId,
                                                    globalPosition.x,
                                                    globalPosition.y,
                                                    globalPosition.z);

      velocity *= 0.0;
      // determine the global coordinates of the next neighbour site:
      // nested loop - x, y, z directions semi-open interval [-2, +2)
      for (site_t x = ((site_t) globalPosition.x) - 1; x < ((site_t) globalPosition.x) + 2; x++)
        for (site_t y = ((site_t) globalPosition.y) - 1; y < ((site_t) globalPosition.y) + 2; y++)
          for (site_t z = ((site_t) globalPosition.z) - 1; z < ((site_t) globalPosition.z) + 2; z++)
          {
            const util::Vector3D<site_t> siteGlobalPosition(x, y, z);

            // convert the global coordinates of the site into a local site index
            proc_t procId;
            site_t siteId;
            bool isValid = latDatLBM->GetContiguousSiteId(siteGlobalPosition, procId, siteId);
            bool isLocal = (procId == topology::NetworkTopology::Instance()->GetLocalRank());

            if (log::Logger::ShouldDisplay<log::Trace>())
            {
              util::Vector3D<site_t> blockCoords, localSiteCoords;
              latDatLBM->GetBlockAndLocalSiteCoords(siteGlobalPosition, blockCoords, localSiteCoords);

              log::Logger::Log<log::Trace, log::OnePerCore>("In colloids::Particle::InterpolateFluidVelocity, particleId: %i,\nsiteGlobalPosition: {%i,%i,%i}\nblockCoords: {%i,%i,%i} - %s\nlocalSiteCoords: {%i,%i,%i} - %s\nprocId: %i, siteId: %i\n",
                                                            particleId,
                                                            siteGlobalPosition.x,
                                                            siteGlobalPosition.y,
                                                            siteGlobalPosition.z,
                                                            blockCoords.x,
                                                            blockCoords.y,
                                                            blockCoords.z,
                                                            !latDatLBM->IsValidBlock(blockCoords)
                                                              ? "INVALID"
                                                              : (latDatLBM->GetBlock(latDatLBM->GetBlockIdFromBlockCoords(blockCoords)).IsEmpty()
                                                                ? "SOLID"
                                                                : "valid"),
                                                            localSiteCoords.x,
                                                            localSiteCoords.y,
                                                            localSiteCoords.z,
                                                            !latDatLBM->IsValidLatticeSite(localSiteCoords)
                                                              ? "INVALID"
                                                              : (!isValid
                                                                ? "SOLID"
                                                                : (!isLocal
                                                                  ? "REMOTE"
                                                                  : "local fluid")),
                                                            procId,
                                                            siteId);
            }

            /** TODO: implement communication for remote fluid sites */
            /** TODO: implement boundary conditions for invalid/solid sites */
            if (!isValid || !isLocal)
              continue;

            // read value of velocity for site index from macroscopic cache
            // TODO: should be LatticeVelocity == Vector3D<LatticeSpeed> (fix as part of #437)
            util::Vector3D<double> siteFluidVelocity = propertyCache.velocityCache.Get(siteId);

            // calculate term of the interpolation sum
            LatticePosition relativePosition(siteGlobalPosition);
            relativePosition -= globalPosition;
            // TODO: should be LatticeVelocity == Vector3D<LatticeSpeed> (fix as part of #437)
            util::Vector3D<double> partialInterpolation = siteFluidVelocity * diracOperation(relativePosition);

            // accumulate each term of the interpolation
            velocity += partialInterpolation;

            log::Logger::Log<log::Trace, log::OnePerCore>("In colloids::Particle::InterpolateFluidVelocity, particleId: %i, siteIndex: %i, fluidVelocity: {%g,%g,%g}, partialInterpolation: {%g,%g,%g}, velocitySoFar: {%g,%g,%g}\n",
                                                          particleId,
                                                          siteId,
                                                          siteFluidVelocity.x,
                                                          siteFluidVelocity.y,
                                                          siteFluidVelocity.z,
                                                          partialInterpolation.x,
                                                          partialInterpolation.y,
                                                          partialInterpolation.z,
                                                          velocity.x,
                                                          velocity.y,
                                                          velocity.z);
          }
      propertyCache.velocityCache.SetRefreshFlag();
    }

  }
}
