// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_KERNELS_GUOFORCINGLBGK_H
#define HEMELB_LB_KERNELS_GUOFORCINGLBGK_H

#include <cstdlib>
#include <cmath>
#include "lb/HFunction.h"
#include "util/utilityFunctions.h"
#include "lb/kernels/BaseKernel.h"
#include "lb/LbmParameters.h"

namespace hemelb::lb
{
    /**
       * Implements the LBGK single-relaxation time kernel, including Guo Forcing
       *
       * Forcing implemented following:
       * Phys. Rev. E 65, 046308 (2002)
       * Zhaoli Guo, Chuguang Zheng, and Baochang Shi
      */
    template<lattice_type LatticeType>
    class GuoForcingLBGK : public BaseKernel<GuoForcingLBGK<LatticeType>, LatticeType>
    {
    public:
        GuoForcingLBGK(InitParams& initParams) :
                BaseKernel<GuoForcingLBGK<LatticeType>, LatticeType>()
        {
        }

        // Adds forcing to momentum
        void DoCalculateDensityMomentumFeq(HydroVars<GuoForcingLBGK>& hydroVars, site_t index);
        // Forwards to LBGK base class
        void DoCalculateFeq(HydroVars<GuoForcingLBGK>&, site_t);
        // Adds forcing to collision term
        void DoCollide(const LbmParameters* const lbmParams,
                       HydroVars<GuoForcingLBGK>& hydroVars);
    };

    template<lattice_type LatticeType>
    struct HydroVars<GuoForcingLBGK<LatticeType> > : HydroVarsBase<LatticeType>
    {
        friend class GuoForcingLBGK<LatticeType>;
    public:
        // Pointer to force at this site
        const LatticeForceVector& force;

        template<class DataSource>
        HydroVars(geometry::Site<DataSource> const &_site) :
                HydroVarsBase<LatticeType>(_site), force(_site.GetForce())
        {
        }

        HydroVars(typename LatticeType::const_span f, const LatticeForceVector& _force) :
                HydroVarsBase<LatticeType>(f), force(_force)
        {
        }

    protected:
        // Guo lattice distribution of external force contributions
        // as calculated in lattice::CalculateForceDistribution.
        FVector<LatticeType> forceDist;
    };

      template<lattice_type LatticeType>
      void GuoForcingLBGK<LatticeType>::DoCalculateDensityMomentumFeq(
          HydroVars<GuoForcingLBGK<LatticeType> >& hydroVars, site_t index)
      {
        LatticeType::CalculateDensityMomentumFEq(hydroVars.f,
                                                 hydroVars.force,
                                                 hydroVars.density,
                                                 hydroVars.momentum,
                                                 hydroVars.velocity,
                                                 hydroVars.f_eq);

        for (unsigned int ii = 0; ii < LatticeType::NUMVECTORS; ++ii)
          hydroVars.f_neq[ii] = hydroVars.f[ii] - hydroVars.f_eq[ii];
      }

      template<lattice_type LatticeType>
      void GuoForcingLBGK<LatticeType>::DoCalculateFeq(HydroVars<GuoForcingLBGK>& hydroVars,
                                                       site_t index)
      {
        LatticeType::CalculateFeq(hydroVars.density,
                                  hydroVars.momentum,
                                  hydroVars.f_eq);

        for (unsigned int ii = 0; ii < LatticeType::NUMVECTORS; ++ii)
          hydroVars.f_neq[ii] = hydroVars.f[ii] - hydroVars.f_eq[ii];
      }

      template<lattice_type LatticeType>
      void GuoForcingLBGK<LatticeType>::DoCollide(const LbmParameters* const lbmParams,
                                                  HydroVars<GuoForcingLBGK>& hydroVars)
      {
        LatticeType::CalculateForceDistribution(lbmParams->GetTau(),
                                                hydroVars.velocity,
                                                hydroVars.force,
                                                hydroVars.forceDist);

        for (Direction dir(0); dir < LatticeType::NUMVECTORS; ++dir)
          hydroVars.SetFPostCollision(dir,
                                      hydroVars.f[dir]
                                          + hydroVars.f_neq[dir] * lbmParams->GetOmega()
                                          + hydroVars.forceDist[dir]);
      }

}

#endif
