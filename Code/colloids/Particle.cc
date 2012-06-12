#include "colloids/Particle.h"

namespace hemelb
{
  namespace colloids
  {
    Particle::Particle(io::xml::XmlAbstractionLayer& xml,
                       lb::MacroscopicPropertyCache& propertyCache) :
      PersistedParticle(xml),
      propertyCache(propertyCache)
    {
    }

    const void Particle::UpdatePosition()
    {
      globalPosition += velocity + bodyForces;
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
            util::Vector3D<int> siteGlobalPosition(x, y, z);

            // convert the global coordinates of the site into a site index
            int siteIndex = 1;//TODO:geometry::ConvertGlobalPositionToSiteIndex(siteGlobalPosition);

            // read value of velocity for site index from macroscopic cache
            util::Vector3D<double> siteFluidVelocity = propertyCache.velocityCache.Get(siteIndex);

            // calculate term of the interpolation sum
            util::Vector3D<double> relativePosition(siteGlobalPosition);
            relativePosition -= globalPosition;
            util::Vector3D<double> partialInterpolation = siteFluidVelocity *
              diracOperation(relativePosition);

            // accumulate each term of the interpolation
            velocity += partialInterpolation;
          }
    }

  }
}
