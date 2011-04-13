#include <limits>
#include "lb/collisions/Collision.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      MinsAndMaxes::MinsAndMaxes()
      {
        MaxDensity = -1.0;
        MaxVelocity = -1.0;
        MaxStress = -1.0;

        MinDensity = std::numeric_limits<double>::max();
        MinVelocity = std::numeric_limits<double>::max();
        MinStress = std::numeric_limits<double>::max();
      }

      // Default constructor, implemented so we can make it protected
      // and prevent instantiation of this base class.
      Collision::Collision()
      {

      }
      Collision::~Collision()
      {
      }

      void Collision::DoCollisions(const bool iDoRayTracing,
                                   const int iFirstIndex,
                                   const int iSiteCount,
                                   const LbmParameters &iLbmParams,
                                   MinsAndMaxes &bMinimaAndMaxima,
                                   geometry::LatticeData &bLatDat,
                                   hemelb::vis::Control *iControl)
      {
        // Standard implementation - do nothing.
      }

      void Collision::PostStep(const bool iDoRayTracing,
                               const int iFirstIndex,
                               const int iSiteCount,
                               const LbmParameters &iLbmParams,
                               MinsAndMaxes &bMinimaAndMaxima,
                               geometry::LatticeData &bLatDat,
                               hemelb::vis::Control *iControl)
      {
        // Standard implementation - do nothing.
      }
    }
  }
}
