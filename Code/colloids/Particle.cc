#include "colloids/Particle.h"

namespace hemelb
{
  namespace colloids
  {
    Particle::Particle(io::xml::XmlAbstractionLayer& xml) :
      PersistedParticle(xml)
    {
    }

    const void Particle::UpdatePosition()
    {
      /** todo: UpdatePosition
       *    - check units
       */
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

    const void Particle::InterpolateFluidVelocity()
    {
      /** todo: InterpolateFluidVelocity
       *    For each local neighbour lattice site
       *    - get velocity for the site from the macroscopic cache object
       *    - calculate site's contribution to the velocity interpolation
       *    - increment particle velocity property with this contribution
       *    - will require communication to transmit remote contributions
       */
      velocity *= 0.0;
    }

  }
}
