#include "colloids/ParticleSet.h"

namespace hemelb
{
  namespace colloids
  {
    ParticleSet::ParticleSet(configuration::XmlAbstractionLayer& xml)
    {
      // assume we are at the <Particles> node
      bool found = xml.MoveToChild("subgridParticle");
      while (found)
      {
        particles.push_back(new Particle(xml));
        found = xml.NextSibling("subgridParticle");
      }
    };
  }
}
