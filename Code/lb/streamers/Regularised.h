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
          typedef typename CollisionType::CKernel::LatticeType LatticeType;

        public:
          Regularised(kernels::InitParams& initParams) :
            collider(initParams)
          {

          }

          template<bool tDoRayTracing>
          inline void DoStreamAndCollide(const site_t firstIndex,
                                         const site_t siteCount,
                                         const LbmParameters* lbmParams,
                                         geometry::LatticeData* latDat,
                                         lb::MacroscopicPropertyCache& propertyCache)
          {
            for (site_t siteIndex = firstIndex; siteIndex < (firstIndex + siteCount); siteIndex++)
            {
              geometry::Site site = latDat->GetSite(siteIndex);

              distribn_t* f = site.GetFOld<LatticeType> ();

              kernels::HydroVars<typename CollisionType::CKernel> hydroVars(f);

              // First calculate the density and macro-velocity
              collider.CalculatePreCollision(hydroVars, site);

              // To evaluate PI, first let unknown particle populations take value given by bounce-back of off-equilibrium parts
              // (fi = fiEq + fopp(i) - fopp(i)Eq)
              distribn_t fTemp[LatticeType::NUMVECTORS];
              for (Direction direction = 0; direction < LatticeType::NUMVECTORS; ++direction)
              {
                if (site.HasBoundary(LatticeType::INVERSEDIRECTIONS[direction]))
                {
                  fTemp[direction] = hydroVars.GetFEq().f[direction] + f[LatticeType::INVERSEDIRECTIONS[direction]]
                      - hydroVars.GetFEq().f[LatticeType::INVERSEDIRECTIONS[direction]];
                }
                else
                {
                  fTemp[direction] = f[direction];
                }
              }

              distribn_t f_neq[LatticeType::NUMVECTORS];
              for (Direction direction = 0; direction < LatticeType::NUMVECTORS; ++direction)
              {
                f_neq[direction] = fTemp[direction] - hydroVars.GetFEq().f[direction];
              }

              // Pi = sum_i e_i e_i f_i
              // zeta = Pi / 2 (Cs^4)
              util::Matrix3D zeta = LatticeType::CalculatePiTensor(f_neq);

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
              site_t streamingDestination[LatticeType::NUMVECTORS];
              for (unsigned l = 0; l < LatticeType::NUMVECTORS; l++)
              {
                if (site.HasBoundary(l))
                {
                  streamingDestination[l] = siteIndex * LatticeType::NUMVECTORS + LatticeType::INVERSEDIRECTIONS[l];
                }
                else
                {
                  streamingDestination[l] = site.GetStreamedIndex<LatticeType> (l);
                }
              }

              const int *Cs[3] = { LatticeType::CX, LatticeType::CY, LatticeType::CZ };

              for (unsigned int ii = 0; ii < LatticeType::NUMVECTORS; ++ii)
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
                    f_neq[ii] += (float (Cs[aa][ii] * Cs[bb][ii])) * zeta[aa][bb];
                  }
                }

                f_neq[ii] *= LatticeType::EQMWEIGHTS[ii];

                /*
                 * Newly constructed distribution function:
                 *    g_i = f^{eq}_i + f^{neq}_i
                 *
                 * Collision step:
                 *    f^{+}_i = g_i + w (g_i - f^{eq}_i)
                 *            = f^{eq}_i + (1+w) f^{neq}_i
                 */
                * (latDat->GetFNew(streamingDestination[ii])) = hydroVars.GetFEq()[ii] + (1.0 + lbmParams->GetOmega())
                    * f_neq[ii];
              }

              ///< @todo #126 It would be nicer if tau is handled in a single place.
              hydroVars.tau = lbmParams->GetTau();

              BaseStreamer<Regularised>::template UpdateMinsAndMaxes<tDoRayTracing>(hydroVars.v_x,
                                                                                    hydroVars.v_y,
                                                                                    hydroVars.v_z,
                                                                                    site,
                                                                                    hydroVars.GetFNeq().f,
                                                                                    hydroVars.density,
                                                                                    hydroVars.tau,
                                                                                    lbmParams,
                                                                                    propertyCache);
            }
          }

          template<bool tDoRayTracing>
          inline void DoPostStep(const site_t iFirstIndex,
                                 const site_t iSiteCount,
                                 const LbmParameters* iLbmParams,
                                 geometry::LatticeData* bLatDat,
                                 lb::MacroscopicPropertyCache& propertyCache)
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
