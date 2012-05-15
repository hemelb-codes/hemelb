#include "colloids/Particle.h"

namespace hemelb
{
  namespace colloids
  {
    Particle::Particle(configuration::XmlAbstractionLayer& xml) :
      PersistedParticle(xml)
    {
    };
  }
}
