#ifndef HEMELB_COLLOIDS_PERSISTEDPARTICLE_H
#define HEMELB_COLLOIDS_PERSISTEDPARTICLE_H

#include "io/xml/XmlAbstractionLayer.h"
#include "util/Vector3D.h"

namespace hemelb
{
  namespace colloids
  {
    class PersistedParticle
    {
      public:
        PersistedParticle(configuration::XmlAbstractionLayer& xml);

      protected:
        PersistedParticle(int particleId, double a0, double ah,
                          util::Vector3D<double> globalPosition) :
          particleId(particleId), smallRadius_a0(a0), largeRadius_ah(ah),
          globalPosition(globalPosition)
        {
        };

        unsigned long          particleId;
        double                 smallRadius_a0;
        double                 largeRadius_ah;
        util::Vector3D<double> globalPosition;
    };
  }
}

#endif /* HEMELB_COLLOIDS_PERSISTEDPARTICLE_H */
