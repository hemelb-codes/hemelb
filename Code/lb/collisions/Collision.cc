#include "lb/collisions/Collision.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      // Default constructor, implemented so we can make it protected
      // and prevent instantiation of this base class.
      Collision::Collision()
      {

      }

      void Collision::DoCollisions(const bool iDoRayTracing,
                                   const double iOmega,
                                   double iFOldAll[],
                                   double iFNewAll[],
                                   const int iFIdAll[],
                                   int iFirstIndex,
                                   const int iSiteCount,
                                   MinsAndMaxes* bMinimaAndMaxima,
                                   const Net* net,
                                   const double iStressType,
                                   const double iStressParam,
                                   hemelb::vis::Control *iControl)
      {
        // Standard implementation - do nothing.
      }

      void Collision::PostStep(const bool iDoRayTracing,
                               const double iOmega,
                               double iFOldAll[],
                               double iFNewAll[],
                               const int iFIdAll[],
                               int iFirstIndex,
                               const int iSiteCount,
                               MinsAndMaxes* bMinimaAndMaxima,
                               const Net* net,
                               const double iStressType,
                               const double iStressParam,
                               hemelb::vis::Control *iControl)
      {
        // Standard implementation - do nothing.
      }
    }
  }
}
