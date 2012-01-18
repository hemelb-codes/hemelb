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
      template<typename CollisionImpl>
      class Regularised : public BaseStreamer<Regularised<CollisionImpl> >
      {
        public:
          typedef CollisionImpl CollisionType;

        private:
          CollisionType collider;

        public:
          Regularised(kernels::InitParams& initParams) :
              collider(initParams)
          {

          }

          template<bool tDoRayTracing>
          inline void DoStreamAndCollide(const site_t iFirstIndex,
                                         const site_t iSiteCount,
                                         const LbmParameters* iLbmParams,
                                         geometry::LatticeData* bLatDat,
                                         hemelb::vis::Control *iControl)
          {
            for (site_t lIndex = iFirstIndex; lIndex < (iFirstIndex + iSiteCount); lIndex++)
            {
              const geometry::Site site = bLatDat->GetSite(lIndex);

              distribn_t* f = site.GetFOld();

              kernels::HydroVars<typename CollisionType::CKernel> hydroVars(f);

              // First calculate the density and macro-velocity
              collider.CalculatePreCollision(hydroVars, site);

              // To evaluate PI, first let unknown particle populations take value given by bounce-back of off-equilibrium parts
              // (fi = fiEq + fopp(i) - fopp(i)Eq)
              distribn_t fTemp[CollisionType::CKernel::LatticeType::NUMVECTORS];
              for (unsigned l = 0; l < CollisionType::CKernel::LatticeType::NUMVECTORS; ++l)
              {
                if (site.HasBoundary(CollisionType::CKernel::LatticeType::INVERSEDIRECTIONS[l]))
                {
                  fTemp[l] = hydroVars.GetFEq().f[l] + f[CollisionType::CKernel::LatticeType::INVERSEDIRECTIONS[l]]
                      - hydroVars.GetFEq().f[CollisionType::CKernel::LatticeType::INVERSEDIRECTIONS[l]];
                }
                else
                {
                  fTemp[l] = f[l];
                }
              }

              distribn_t f_neq[CollisionType::CKernel::LatticeType::NUMVECTORS];
              for (unsigned l = 0; l < CollisionType::CKernel::LatticeType::NUMVECTORS; ++l)
              {
                f_neq[l] = fTemp[l] - hydroVars.GetFEq().f[l];
              }

              // Pi = sum_i e_i e_i f_i
              // zeta = Pi / 2 (Cs^4)
              Order2Tensor zeta = CollisionType::CKernel::LatticeType::CalculatePiTensor(f_neq);

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
              site_t lStreamTo[CollisionType::CKernel::LatticeType::NUMVECTORS];
              for (unsigned l = 0; l < CollisionType::CKernel::LatticeType::NUMVECTORS; l++)
              {
                if (site.HasBoundary(l))
                {
                  lStreamTo[l] = lIndex * CollisionType::CKernel::LatticeType::NUMVECTORS
                      + CollisionType::CKernel::LatticeType::INVERSEDIRECTIONS[l];
                }
                else
                {
                  lStreamTo[l] = site.GetStreamedIndex(l);
                }
              }

              const int *Cs[3] = { CollisionType::CKernel::LatticeType::CX, CollisionType::CKernel::LatticeType::CY,
                                   CollisionType::CKernel::LatticeType::CZ };

              for (unsigned int ii = 0; ii < CollisionType::CKernel::LatticeType::NUMVECTORS; ++ii)
              {
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

                f_neq[ii] *= CollisionType::CKernel::LatticeType::EQMWEIGHTS[ii];

                /*
                 * Newly constructed distribution function:
                 *    g_i = f^{eq}_i + f^{neq}_i
                 *
                 * Collision step:
                 *    f^{+}_i = g_i + w (g_i - f^{eq}_i)
                 *            = f^{eq}_i + (1+w) f^{neq}_i
                 */
                * (bLatDat->GetFNew(lStreamTo[ii])) = hydroVars.GetFEq().f[ii]
                    + (1.0 + iLbmParams->GetOmega()) * f_neq[ii];
              }

              BaseStreamer<Regularised>::template UpdateMinsAndMaxes<tDoRayTracing>(hydroVars.v_x,
                                                                                    hydroVars.v_y,
                                                                                    hydroVars.v_z,
                                                                                    site,
                                                                                    hydroVars.GetFNeq().f,
                                                                                    hydroVars.density,
                                                                                    iLbmParams,
                                                                                    iControl);
            }
          }

          template<bool tDoRayTracing>
          inline void DoPostStep(const site_t iFirstIndex,
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
