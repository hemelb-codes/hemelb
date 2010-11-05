#ifndef HEMELB_LB_COLLISIONS_COLLISION_H
#define HEMELB_LB_COLLISIONS_COLLISION_H

#include "vis/Control.h"
#include "lb/LocalLatticeData.h"

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
                                    const int iFirstIndex,
                                    const int iSiteCount,
                                    MinsAndMaxes* bMinimaAndMaxima,
                                    LocalLatticeData &bLocalLatDat,
                                    const double iStressType,
                                    const double iStressParam,
                                    hemelb::vis::Control *iControl);

          virtual void PostStep(const bool iDoRayTracing,
                                const double iOmega,
                                const int iFirstIndex,
                                const int iSiteCount,
                                MinsAndMaxes* bMinimaAndMaxima,
                                LocalLatticeData &bLocalLatDat,
                                const double iStressType,
                                const double iStressParam,
                                hemelb::vis::Control *iControl);

        protected:
          // Use a protected constructor to ensure the class is never instantiated.
          Collision();
          virtual ~Collision();

          template<bool tDoRayTracing>
          void UpdateMinsAndMaxes(double iVx,
                                  double iVy,
                                  double iVz,
                                  const int &iSiteIndex,
                                  const double *f_neq,
                                  const double &iDensity,
                                  MinsAndMaxes *bMinimaAndMaxima,
                                  const LocalLatticeData &iLocalLatDat,
                                  const double &iStressType,
                                  const double &iStressParam,
                                  hemelb::vis::Control *iControl)
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
              if (iLocalLatDat.GetNormalToWall(iSiteIndex)[0] > 1.0e+30)
              {
                stress = 0.0;
                rtStress = 1.0e+30;
              }
              else
              {
                D3Q15::CalculateShearStress(
                                            iDensity,
                                            f_neq,
                                            iLocalLatDat.GetNormalToWall(
                                                                         iSiteIndex),
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
              iControl->RegisterSite(iSiteIndex, iDensity, lVelocity, rtStress);
            }
          }
      };
    }
  }
}

#endif //HEMELB_LB_COLLISIONS_COLLISION_H
