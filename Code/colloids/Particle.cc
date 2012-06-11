#include "colloids/Particle.h"

namespace hemelb
{
  namespace colloids
  {
    Particle::Particle(io::xml::XmlAbstractionLayer& xml) :
      PersistedParticle(xml)
    {
    };
  }
}
