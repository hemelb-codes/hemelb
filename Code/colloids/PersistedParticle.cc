#include "colloids/PersistedParticle.h"

namespace hemelb
{
  namespace colloids
  {
    PersistedParticle::PersistedParticle(io::xml::XmlAbstractionLayer& xml)
    {
      // assume we are currently at a <SubgridParticle> node

      bool ok = true;
      ok &= xml.GetUnsignedLongValue("ParticleId", particleId);
      ok &= xml.GetDoubleValue("InputRadiusA0", smallRadius_a0);
      ok &= xml.GetDoubleValue("HydrostaticRadiusAh", largeRadius_ah);
      ok &= xml.GetDoubleVector("initialPosition", globalPosition);
    };
  }
}
