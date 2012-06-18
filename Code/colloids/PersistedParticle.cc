#include "colloids/PersistedParticle.h"
#include "log/Logger.h"

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

      /** TODO: convert units from physical distance to lattice distance */

      if (log::Logger::ShouldDisplay<log::Debug>())
        log::Logger::Log<log::Debug, log::OnePerCore>(
          "In colloids::PersistedParticle::ctor, id: %i, a0: %g, ah: %g, position: {%g,%g,%g}\n",
          particleId, smallRadius_a0, largeRadius_ah,
          globalPosition.x, globalPosition.y, globalPosition.z);
    };

  }
}
