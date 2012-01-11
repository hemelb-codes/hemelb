#ifndef HEMELB_LB_STREAMERS_REGULARISED_H
#define HEMELB_LB_STREAMERS_REGULARISED_H

#include "lb/streamers/BaseStreamer.h"

namespace hemelb
{
  namespace lb
  {
    namespace streamers
    {
      // TODO REFACTOR this class to be just a collision, using the BounceBack streamer.
      template<typename CollisionType>
      class Regularised : public BaseStreamer<Regularised<CollisionType> >
      {
        private:
          CollisionType collider;

        public:
          Regularised(kernels::InitParams& initParams) :
              collider(initParams)
          {

          }

          template<bool tDoRayTracing>
          void DoStreamAndCollide(const site_t iFirstIndex,
                                  const site_t iSiteCount,
                                  const LbmParameters* iLbmParams,
                                  geometry::LatticeData* bLatDat,
                                  hemelb::vis::Control *iControl)
          {
            for (site_t lIndex = iFirstIndex; lIndex < (iFirstIndex + iSiteCount); lIndex++)
            {
              distribn_t* f = bLatDat->GetFOld(lIndex * D3Q15::NUMVECTORS);

              kernels::HydroVars<typename CollisionType::CKernel> hydroVars(f);

              // First calculate the density and macro-velocity
              //! @todo does the following comment make sense anymore?
              // TEMPORARILY STORE f_eq IN f_neq BUT THE FUNCTION RETURNS f_eq. THIS IS SORTED
              // OUT IN A SUBSEQUENT FOR LOOP.
              collider.CalculatePreCollision(hydroVars, lIndex - iFirstIndex);

              // To evaluate PI, first let unknown particle populations take value given by bounce-back of off-equilibrium parts
              // (fi = fiEq + fopp(i) - fopp(i)Eq)
              distribn_t fTemp[D3Q15::NUMVECTORS];

              for (unsigned l = 0; l < D3Q15::NUMVECTORS; ++l)
              {
                //! @todo I thought the components that need updating are those that satisfy bLatDat->HasBoundary(lIndex, D3Q15::INVERSEDIRECTIONS[l]), i.e. not having a node streaming on them.
                if (bLatDat->HasBoundary(lIndex, l))
                {
                  //! @todo please document how fi = fiEq + fopp(i) - fopp(i)Eq becomes the expression below
                  fTemp[l] = f[D3Q15::INVERSEDIRECTIONS[l]]
                      + 3.0 * D3Q15::EQMWEIGHTS[l]
                          * (hydroVars.v_x * D3Q15::CX[l] + hydroVars.v_y * D3Q15::CY[l]
                              + hydroVars.v_z * D3Q15::CZ[l]);
                }
                else
                {
                  fTemp[l] = f[l];
                }
              }

              distribn_t f_neq[D3Q15::NUMVECTORS];
              for (unsigned l = 0; l < D3Q15::NUMVECTORS; ++l)
              {
                f_neq[l] = fTemp[l] - hydroVars.GetFEq().f[l];
              }

              // Pi = sum_i e_i e_i f_i
              // zeta = Pi / 2 (Cs^4)
              Order2Tensor zeta = D3Q15::CalculatePiTensor(f_neq);

              for (int m = 0; m < 3; m++)
              {
                for (int n = 0; n < 3; n++)
                {
                  zeta[m][n] /= (2.0 * Cs2 * Cs2);
                }
              }

              // chi = Cs^2 I : zeta
              const distribn_t chi = Cs2 * (zeta[0][0] + zeta[1][1] + zeta[2][2]);

              // Now apply bounce-back to the components that require it
              site_t lStreamTo[D3Q15::NUMVECTORS];
              for (unsigned l = 0; l < D3Q15::NUMVECTORS; l++)
              {
                if (bLatDat->HasBoundary(lIndex, l))
                {
                  lStreamTo[l] = lIndex * D3Q15::NUMVECTORS + D3Q15::INVERSEDIRECTIONS[l];
                }
                else
                {
                  lStreamTo[l] = bLatDat->GetStreamedIndex(lIndex, l);
                }
              }

              const int *Cs[3] = { D3Q15::CX, D3Q15::CY, D3Q15::CZ };

              for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ++ii)
              {
                // Calculate the dot-product of the velocity with the direction vector.
                distribn_t vSum = hydroVars.v_x * (float) D3Q15::CX[ii]
                    + hydroVars.v_y * (float) D3Q15::CY[ii] + hydroVars.v_z * (float) D3Q15::CZ[ii];

                // Calculate the squared magnitude of the velocity.
                distribn_t v2Sum = hydroVars.v_x * hydroVars.v_x + hydroVars.v_y * hydroVars.v_y
                    + hydroVars.v_z * hydroVars.v_z;

                // F eqm = density proportional component...
                distribn_t streamed = hydroVars.density;

                // ... - v^2 component...
                streamed -= ( (3.0 / 2.0) * v2Sum / hydroVars.density);

                // ... + v^1 component
                streamed += 3.0 * vSum + (9.0 / 2.0) * vSum * vSum / hydroVars.density;

                // Multiply by eqm weight.
                streamed *= D3Q15::EQMWEIGHTS[ii];

                //! @todo: I think "streamed" is already known as hydroVars.GetFEq().f[ii]

                // According to Latt & Chopard (Physical Review E77, 2008),
                // f_neq[i] = (LatticeWeight[i] / (2 Cs^4)) *
                //            Q_i : Pi(n_eq)
                // Where Q_i = c_i c_i - Cs^2 I
                // and Pi(n_eq) = Sum{i} (c_i c_i f_i)
                //
                // We pre-compute zeta = Pi(neq) / (2 Cs^4)
                //             and chi =  Cs^2 I : zeta
                // Hence we can compute f_neq[i] = LatticeWeight[i] * ((c_i c_i) : zeta - chi)
                f_neq[ii] = -chi;

                for (int aa = 0; aa < 3; ++aa)
                {
                  for (int bb = 0; bb < 3; ++bb)
                  {
                    f_neq[ii] += (float(Cs[aa][ii] * Cs[bb][ii])) * zeta[aa][bb];
                  }
                }

                f_neq[ii] *= D3Q15::EQMWEIGHTS[ii];

                /*
                 * Newly constructed distribution function:
                 *    g_i = f^{eq}_i + f^{neq}_i
                 *
                 * Collision step:
                 *    f^{+}_i = g_i + w (g_i - f^{eq}_i)
                 *            = f^{eq}_i + (1+w) f^{neq}_i
                 */
                * (bLatDat->GetFNew(lStreamTo[ii])) = streamed
                    + (1.0 + iLbmParams->GetOmega()) * f_neq[ii];
              }

              BaseStreamer<Regularised>::template UpdateMinsAndMaxes<tDoRayTracing>(hydroVars.v_x,
                                                                                    hydroVars.v_y,
                                                                                    hydroVars.v_z,
                                                                                    lIndex,
                                                                                    hydroVars.GetFNeq().f,
                                                                                    hydroVars.density,
                                                                                    bLatDat,
                                                                                    iLbmParams,
                                                                                    iControl);
            }
          }

          template<bool tDoRayTracing>
          void DoPostStep(const site_t iFirstIndex,
                          const site_t iSiteCount,
                          const LbmParameters* iLbmParams,
                          geometry::LatticeData* bLatDat,
                          hemelb::vis::Control *iControl)
          {

          }

          void DoReset(kernels::InitParams* init)
          {
            collider.Reset(init);
          }
      };
    }
  }
}

#endif /* HEMELB_LB_STREAMERS_REGULARISED_H */
