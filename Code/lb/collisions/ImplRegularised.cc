#include "lb/collisions/ImplRegularised.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      void ImplRegularised::DoCollisions(const bool iDoRayTracing,
                                         const int iFirstIndex,
                                         const int iSiteCount,
                                         const LbmParameters &iLbmParams,
                                         MinsAndMaxes &bMinimaAndMaxima,
                                         geometry::LatticeData &bLatDat,
                                         hemelb::vis::Control *iControl)
      {
        if (iDoRayTracing)
        {
          DoCollisionsInternal<true> (iFirstIndex, iSiteCount, iLbmParams, bMinimaAndMaxima,
                                      bLatDat, iControl);
        }
        else
        {
          DoCollisionsInternal<false> (iFirstIndex, iSiteCount, iLbmParams, bMinimaAndMaxima,
                                       bLatDat, iControl);
        }
      }

      template<bool tDoRayTracing>
      void ImplRegularised::DoCollisionsInternal(const int iFirstIndex,
                                                 const int iSiteCount,
                                                 const LbmParameters &iLbmParams,
                                                 MinsAndMaxes &bMinimaAndMaxima,
                                                 geometry::LatticeData &bLatDat,
                                                 hemelb::vis::Control *iControl)
      {
        for (int lIndex = iFirstIndex; lIndex < (iFirstIndex + iSiteCount); lIndex++)
        {
          double *f = bLatDat.GetFOld(lIndex * D3Q15::NUMVECTORS);
          double density, v_x, v_y, v_z;
          double f_neq[15];

          // First calculate the density and macro-velocity
          // TEMPORARILY STORE f_eq IN f_neq BUT THE FUNCTION RETURNS f_eq. THIS IS SORTED
          // OUT IN A SUBSEQUENT FOR LOOP.
          D3Q15::CalculateDensityVelocityFEq(f, density, v_x, v_y, v_z, f_neq);

          // To evaluate PI, first let unknown particle populations take value given by bounce-back of off-equilibrium parts
          // (fi = fiEq + fopp(i) - fopp(i)Eq)
          double fTemp[15];

          for (int l = 0; l < 15; l++)
            fTemp[l] = f[l];

          if (bLatDat.HasBoundary(lIndex, 1))
          {
            fTemp[1] = f[2] + (2.0 / 3.0) * v_x;
          }
          if (bLatDat.HasBoundary(lIndex, 2))
          {
            fTemp[2] = f[1] - (2.0 / 3.0) * v_x;
          }
          if (bLatDat.HasBoundary(lIndex, 3))
          {
            fTemp[3] = f[4] + (2.0 / 3.0) * v_y;
          }
          if (bLatDat.HasBoundary(lIndex, 4))
          {
            fTemp[4] = f[3] - (2.0 / 3.0) * v_y;
          }
          if (bLatDat.HasBoundary(lIndex, 5))
          {
            fTemp[5] = f[6] + (2.0 / 3.0) * v_z;
          }
          if (bLatDat.HasBoundary(lIndex, 6))
          {
            fTemp[6] = f[5] - (2.0 / 3.0) * v_z;
          }
          if (bLatDat.HasBoundary(lIndex, 7))
          {
            fTemp[7] = f[8] + (2.0 / 24.0) * ( (v_x + v_y) + v_z);
          }
          if (bLatDat.HasBoundary(lIndex, 8))
          {
            fTemp[8] = f[7] - (2.0 / 24.0) * ( (v_x + v_y) + v_z);
          }
          if (bLatDat.HasBoundary(lIndex, 9))
          {
            fTemp[9] = f[10] + (2.0 / 24.0) * ( (v_x + v_y) - v_z);
          }
          if (bLatDat.HasBoundary(lIndex, 10))
          {
            fTemp[10] = f[9] - (2.0 / 24.0) * ( (v_x + v_y) - v_z);
          }
          if (bLatDat.HasBoundary(lIndex, 11))
          {
            fTemp[11] = f[12] + (2.0 / 24.0) * ( (v_x - v_y) + v_z);
          }
          if (bLatDat.HasBoundary(lIndex, 12))
          {
            fTemp[12] = f[11] - (2.0 / 24.0) * ( (v_x - v_y) + v_z);
          }
          if (bLatDat.HasBoundary(lIndex, 13))
          {
            fTemp[13] = f[14] + (2.0 / 24.0) * ( (v_x - v_y) - v_z);
          }
          if (bLatDat.HasBoundary(lIndex, 14))
          {
            fTemp[14] = f[13] - (2.0 / 24.0) * ( (v_x - v_y) - v_z);
          }

          // UP TO THIS POINT, F_NEQ ACTUALLY CONTAINS F_EQ. AT THIS
          // STAGE WE REPLACE IT WITH THE ACTUAL NON-EQ VALUE, POST
          // BOUNCING_BACK WHERE NEEDED.
          for (int l = 0; l < 15; l++)
          {
            f_neq[l] = fTemp[l] - f_neq[l];
          }

          double density_1 = 1. / density;
          double v_xx = v_x * v_x;
          double v_yy = v_y * v_y;
          double v_zz = v_z * v_z;

          // PI = sum_i e_i e_i f_i
          double piMatrix[3][3];

          double diagSum = f_neq[7] + f_neq[8] + f_neq[9] + f_neq[10] + f_neq[11] + f_neq[12]
              + f_neq[13] + f_neq[14];

          piMatrix[0][0] = f_neq[1] + f_neq[2] + diagSum;
          piMatrix[0][1] = f_neq[7] + f_neq[8] + f_neq[9] + f_neq[10] - (f_neq[11] + f_neq[12]
              + f_neq[13] + f_neq[14]);
          piMatrix[0][2] = f_neq[7] + f_neq[8] + f_neq[11] + f_neq[12] - (f_neq[9] + f_neq[10]
              + f_neq[13] + f_neq[14]);
          piMatrix[1][0] = piMatrix[0][1];
          piMatrix[1][1] = f_neq[3] + f_neq[4] + diagSum;
          piMatrix[1][2] = f_neq[7] + f_neq[8] + f_neq[13] + f_neq[14] - (f_neq[9] + f_neq[10]
              + f_neq[11] + f_neq[12]);
          piMatrix[2][0] = piMatrix[0][2];
          piMatrix[2][1] = piMatrix[1][2];
          piMatrix[2][2] = f_neq[5] + f_neq[6] + diagSum;

          for (int m = 0; m < 3; m++)
            for (int n = 0; n < 3; n++)
              piMatrix[m][n] /= (2.0 * Cs2 * Cs2);

          // Qi = e_i e_i - (speed of sound ^ 2) * Identity
          // Then gi = fiEq + t_i (the 2/9, 1/9, 1/72 stuff) (Qi . PI (inner product)) / 2 * speed of sound^4
          // Or:  gi = fiEq + t_i (the 2/9, 1/9, 1/72 stuff) ((e_i e_i . PI (inner product)) / 2 * speed of sound^4 - specialNumber)
          double specialNumber = (2.0 / 9.0) * Cs2 * (piMatrix[0][0] + piMatrix[1][1]
              + piMatrix[2][2]);
          double piMatrixSum = piMatrix[0][0] + piMatrix[0][1] + piMatrix[0][2] + piMatrix[1][0]
              + piMatrix[1][1] + piMatrix[1][2] + piMatrix[2][0] + piMatrix[2][1] + piMatrix[2][2];

          // The gi (here; f) are then collided and streamed
          * (bLatDat.GetFNew(bLatDat.GetStreamedIndex(lIndex, 0))) = ( (2.0 / 9.0) * density - (1.0
              / 3.0) * ( (v_xx + v_yy + v_zz) * density_1)) + (1.0 + iLbmParams.Omega) * (f_neq[0]
              = -specialNumber);

          double temp1 = (1.0 / 9.0) * density - (1.0 / 6.0) * ( (v_xx + v_yy + v_zz) * density_1);
          specialNumber *= 1.0 / 2.0;

          // Now apply bounce-back to the components that require it, from fTemp
          int lStreamTo[15];
          for (int l = 1; l < 15; l++)
          {
            if (bLatDat.HasBoundary(lIndex, l))
            {
              lStreamTo[l] = lIndex * 15 + D3Q15::INVERSEDIRECTIONS[l];
            }
            else
            {
              lStreamTo[l] = bLatDat.GetStreamedIndex(lIndex, l);
            }
          }

          * (bLatDat.GetFNew(lStreamTo[1])) = temp1 + (0.5 * density_1) * v_xx + (1.0 / 3.0) * v_x
              + (1.0 + iLbmParams.Omega)
                  * (f_neq[1] = (1.0 / 9.0) * piMatrix[0][0] - specialNumber); // (+1, 0, 0)
          * (bLatDat.GetFNew(lStreamTo[2])) = temp1 + (0.5 * density_1) * v_xx - (1.0 / 3.0) * v_x
              + (1.0 + iLbmParams.Omega)
                  * (f_neq[2] = (1.0 / 9.0) * piMatrix[0][0] - specialNumber); // (+1, 0, 0)

          * (bLatDat.GetFNew(lStreamTo[3])) = temp1 + (0.5 * density_1) * v_yy + (1.0 / 3.0) * v_y
              + (1.0 + iLbmParams.Omega)
                  * (f_neq[3] = (1.0 / 9.0) * piMatrix[1][1] - specialNumber); // (0, +1, 0)
          * (bLatDat.GetFNew(lStreamTo[4])) = temp1 + (0.5 * density_1) * v_yy - (1.0 / 3.0) * v_y
              + (1.0 + iLbmParams.Omega)
                  * (f_neq[4] = (1.0 / 9.0) * piMatrix[1][1] - specialNumber); // (0, +1, 0)

          * (bLatDat.GetFNew(lStreamTo[5])) = temp1 + (0.5 * density_1) * v_zz + (1.0 / 3.0) * v_z
              + (1.0 + iLbmParams.Omega)
                  * (f_neq[5] = (1.0 / 9.0) * piMatrix[2][2] - specialNumber); // (0, +1, 0)
          * (bLatDat.GetFNew(lStreamTo[6])) = temp1 + (0.5 * density_1) * v_zz - (1.0 / 3.0) * v_z
              + (1.0 + iLbmParams.Omega)
                  * (f_neq[6] = (1.0 / 9.0) * piMatrix[2][2] - specialNumber); // (0, +1, 0)

          temp1 *= (1.0 / 8.0);
          specialNumber *= (1.0 / 8.0);

          double temp2 = (v_x + v_y) + v_z;

          * (bLatDat.GetFNew(lStreamTo[7])) = temp1 + (1.0 / 16.0) * density_1 * temp2 * temp2
              + (1.0 / 24.0) * temp2 + (1.0 + iLbmParams.Omega) * (f_neq[7] = ( (1.0 / 72.0)
              * piMatrixSum - specialNumber)); // (+1, +1, +1)
          * (bLatDat.GetFNew(lStreamTo[8])) = temp1 + (1.0 / 16.0) * density_1 * temp2 * temp2
              + (-1.0 / 24.0) * temp2 + (1.0 + iLbmParams.Omega) * (f_neq[8] = ( (1.0 / 72.0)
              * piMatrixSum - specialNumber)); // (-1, -1, -1)

          temp2 = (v_x + v_y) - v_z;

          * (bLatDat.GetFNew(lStreamTo[9])) = temp1 + (1.0 / 16.0) * density_1 * temp2 * temp2
              + (1.0 / 24.0) * temp2 + (1.0 + iLbmParams.Omega) * (f_neq[9] = ( (1.0 / 72.0)
              * (piMatrixSum - 4.0 * (piMatrix[0][2] + piMatrix[1][2])) - specialNumber)); // (+1, +1, -1)
          * (bLatDat.GetFNew(lStreamTo[10])) = temp1 + (1.0 / 16.0) * density_1 * temp2 * temp2
              + (-1.0 / 24.0) * temp2 + (1.0 + iLbmParams.Omega) * (f_neq[10] = ( (1.0 / 72.0)
              * (piMatrixSum - 4.0 * (piMatrix[0][2] + piMatrix[1][2])) - specialNumber)); // (-1, -1, +1)

          temp2 = (v_x - v_y) + v_z;

          * (bLatDat.GetFNew(lStreamTo[11])) = temp1 + (1.0 / 16.0) * density_1 * temp2 * temp2
              + (1.0 / 24.0) * temp2 + (1.0 + iLbmParams.Omega) * (f_neq[11] = ( (1.0 / 72.0)
              * (piMatrixSum - 4.0 * (piMatrix[0][1] + piMatrix[1][2])) - specialNumber)); // (+1, -1, +1)
          * (bLatDat.GetFNew(lStreamTo[12])) = temp1 + (1.0 / 16.0) * density_1 * temp2 * temp2
              + (-1.0 / 24.0) * temp2 + (1.0 + iLbmParams.Omega) * (f_neq[12] = ( (1.0 / 72.0)
              * (piMatrixSum - 4.0 * (piMatrix[0][1] + piMatrix[1][2])) - specialNumber)); // (-1, +1, -1)

          temp2 = (v_x - v_y) - v_z;

          * (bLatDat.GetFNew(lStreamTo[13])) = temp1 + (1.0 / 16.0) * density_1 * temp2 * temp2
              + (1.0 / 24.0) * temp2 + (1.0 + iLbmParams.Omega) * (f_neq[13] = ( (1.0 / 72.0)
              * (piMatrixSum - 4.0 * (piMatrix[0][1] + piMatrix[0][2])) - specialNumber)); // (+1, -1, -1)
          * (bLatDat.GetFNew(lStreamTo[14])) = temp1 + (1.0 / 16.0) * density_1 * temp2 * temp2
              + (-1.0 / 24.0) * temp2 + (1.0 + iLbmParams.Omega) * (f_neq[14] = ( (1.0 / 72.0)
              * (piMatrixSum - 4.0 * (piMatrix[0][1] + piMatrix[0][2])) - specialNumber)); // (-1, +1, +1)

          UpdateMinsAndMaxes<tDoRayTracing> (v_x, v_y, v_z, lIndex, f_neq, density,
                                             bMinimaAndMaxima, bLatDat, iLbmParams, iControl);
        }
      }
    }
  }
}
