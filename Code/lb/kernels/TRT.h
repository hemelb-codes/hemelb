
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_KERNELS_TRT_H
#define HEMELB_LB_KERNELS_TRT_H

#include <cstdlib>
#include "lb/HFunction.h"
#include "util/utilityFunctions.h"
#include "lb/kernels/BaseKernel.h"

namespace hemelb
{
  namespace lb
  {
    namespace kernels
    {
      /**
       * TRT: This class implements a two-relaxation time kernel.
       */
      template<class LatticeType>
      class TRT : public BaseKernel<TRT<LatticeType>, LatticeType>
      {
          // Store the directions as pairs of opposites
          // (Note that the zero vector is it's own opposite)
          typedef std::pair<Direction, Direction> Opposites;
          typedef std::vector<Opposites> OppList;
          OppList directionPairs;
          Direction iZero;

        public:
          TRT(InitParams& initParams)
          {
            for (Direction i = 0; i < LatticeType::NUMVECTORS; ++i)
            {
              Direction iBar = LatticeType::INVERSEDIRECTIONS[i];
              if (i == iBar)
                iZero = i;

              if (iBar >= i)
                directionPairs.push_back(std::make_pair(i, iBar));
            }
          }

          inline void DoCalculateDensityMomentumFeq(HydroVars<TRT<LatticeType> >& hydroVars, site_t index)
          {
            LatticeType::CalculateDensityMomentumFEq(hydroVars.f,
                                                     hydroVars.density,
                                                     hydroVars.momentum.x,
                                                     hydroVars.momentum.y,
                                                     hydroVars.momentum.z,
                                                     hydroVars.velocity.x,
                                                     hydroVars.velocity.y,
                                                     hydroVars.velocity.z,
                                                     hydroVars.f_eq.f);

            for (unsigned int ii = 0; ii < LatticeType::NUMVECTORS; ++ii)
            {
              hydroVars.f_neq.f[ii] = hydroVars.f[ii] - hydroVars.f_eq.f[ii];
            }
          }

          inline void DoCalculateFeq(HydroVars<TRT>& hydroVars, site_t index)
          {
            LatticeType::CalculateFeq(hydroVars.density,
                                      hydroVars.momentum.x,
                                      hydroVars.momentum.y,
                                      hydroVars.momentum.z,
                                      hydroVars.f_eq.f);

            for (unsigned int ii = 0; ii < LatticeType::NUMVECTORS; ++ii)
            {
              hydroVars.f_neq.f[ii] = hydroVars.f[ii] - hydroVars.f_eq.f[ii];
            }
          }

          inline void DoCollide(const LbmParameters* const lbmParams, HydroVars<TRT>& hydroVars)
          {
            // Note HemeLB defines omega = -1/ tau
            // Magic number determines the other relaxation time
            // Lambda = (tau_plus - 1/2) (tau_minus - 1/2)
            // Choose such that HWBB walls are always in the right place.
            // TODO: make this a configurable parameter.
            const distribn_t Lambda = 3.0 / 16.0;

            const distribn_t tau_plus = lbmParams->GetTau();
            const distribn_t omega_plus = lbmParams->GetOmega();
            const distribn_t tau_minus = 0.5 + Lambda / (tau_plus - 0.5);
            const distribn_t omega_minus =  -1.0 / tau_minus;

            // Special case the null velocity.
            hydroVars.SetFPostCollision(iZero,
                                        hydroVars.f[iZero] + omega_plus * hydroVars.f_neq.f[iZero]);

            // Now deal with the non-zero
            for (OppList::const_iterator oppIt = directionPairs.begin();
                oppIt != directionPairs.end();
                ++oppIt)
            {
              const Direction& i = oppIt->first;
              const Direction& iBar = oppIt->second;

              distribn_t sym = 0.5 * omega_plus * (hydroVars.f_neq.f[i] + hydroVars.f_neq.f[iBar]);
              distribn_t asym = 0.5 * omega_minus * (hydroVars.f_neq.f[i] - hydroVars.f_neq.f[iBar]);
              hydroVars.SetFPostCollision(i, hydroVars.f[i] + sym + asym);
              hydroVars.SetFPostCollision(iBar, hydroVars.f[iBar] + sym - asym);
            }
          }

      };

    }
  }
}

#endif /* HEMELB_LB_KERNELS_TRT_H */
