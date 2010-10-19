#ifndef HEMELB_LB_COLLISIONS_COLLISION_H
#define HEMELB_LB_COLLISIONS_COLLISION_H

// TODO BOO AND HISS. Find some other way.
#include "net.h"
#include "vis/RayTracer.h"

#include <math.h>

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      struct MinsAndMaxes
      {
        public:
          double MinDensity;
          double MaxDensity;
          double MinVelocity;
          double MaxVelocity;
          double MinStress;
          double MaxStress;
      };

      class Collision
      {
        public:
          virtual void DoCollisions(const bool iDoRayTracing,
                                    const double iOmega,
                                    double iFOldAll[],
                                    double iFNewAll[],
                                    const int iFIdAll[],
                                    const int iFirstIndex,
                                    const int iSiteCount,
                                    MinsAndMaxes* bMinimaAndMaxima,
                                    const Net* net,
                                    const double iStressType,
                                    const double iStressParam);

          virtual void PostStep(const bool iDoRayTracing,
                                const double iOmega,
                                double iFOldAll[],
                                double iFNewAll[],
                                const int iFIdAll[],
                                const int iFirstIndex,
                                const int iSiteCount,
                                MinsAndMaxes* bMinimaAndMaxima,
                                const Net* net,
                                const double iStressType,
                                const double iStressParam);

        protected:
          // Use a protected constructor to ensure the class is never instantiated.
          Collision();

          template<bool tDoRayTracing>
          void UpdateMinsAndMaxes(double iVx,
                                  double iVy,
                                  double iVz,
                                  const int &iSiteIndex,
                                  const double *f_neq,
                                  const double &iDensity,
                                  MinsAndMaxes *bMinimaAndMaxima,
                                  const Net* net,
                                  const double &iStressType,
                                  const double &iStressParam)
          {
            if (iDensity < bMinimaAndMaxima->MinDensity)
            {
              bMinimaAndMaxima->MinDensity = iDensity;
            }
            if (iDensity > bMinimaAndMaxima->MaxDensity)
            {
              bMinimaAndMaxima->MaxDensity = iDensity;
            }

            double stress;
            double rtStress;

            // TODO: It'd be nice if this were unnecessary.
            iVx *= (1.0 / iDensity);
            iVy *= (1.0 / iDensity);
            iVz *= (1.0 / iDensity);

            double lVelocity = sqrt(iVx * iVx + iVy * iVy + iVz * iVz);

            if (lVelocity < bMinimaAndMaxima->MinVelocity)
            {
              bMinimaAndMaxima->MinVelocity = lVelocity;
            }
            if (lVelocity > bMinimaAndMaxima->MaxVelocity)
            {
              bMinimaAndMaxima->MaxVelocity = lVelocity;
            }

            if (iStressType == SHEAR_STRESS)
            {
              if (net->GetNormalToWall(iSiteIndex)[0] > 1.0e+30)
              {
                stress = 0.0;
                rtStress = 1.0e+30;
              }
              else
              {
                D3Q15::CalculateShearStress(iDensity, f_neq,
                                            net->GetNormalToWall(iSiteIndex),
                                            stress, iStressParam);
                rtStress = stress;
              }
            }
            else
            {
              D3Q15::CalculateVonMisesStress(f_neq, stress, iStressParam);
              rtStress = stress;
            }

            if (stress < bMinimaAndMaxima->MinStress)
            {
              bMinimaAndMaxima->MinStress = stress;
            }
            if (stress > bMinimaAndMaxima->MaxStress)
            {
              bMinimaAndMaxima->MaxStress = stress;
            }

            if (tDoRayTracing)
            {
              hemelb::vis::rtUpdateClusterVoxel(iSiteIndex, iDensity,
                                                lVelocity, rtStress);
            }
          }
      };
    }
  }
}

#endif //HEMELB_LB_COLLISIONS_COLLISION_H
