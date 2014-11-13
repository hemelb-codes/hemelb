// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_LB_KERNELS_GUO_FORCING_LBGK_H
#define HEMELB_LB_KERNELS_GUO_FORCING_LBGK_H

#include <cstdlib>
#include <cmath>
#include "lb/HFunction.h"
#include "util/utilityFunctions.h"
#include "lb/kernels/BaseKernel.h"

namespace hemelb { namespace lb { namespace kernels {
  /**
   * Implements the LBGK single-relaxation time kernel, including Guo Forcing
   *
   * Forcing implemented following:
   * Phys. Rev. E 65, 046308 (2002)
   * Zhaoli Guo, Chuguang Zheng, and Baochang Shi
   */
  template<class LatticeType>
    class GuoForcingLBGK
        : public BaseKernel<GuoForcingLBGK<LatticeType>, LatticeType> {
      public:
        GuoForcingLBGK(InitParams& initParams)
          : BaseKernel<GuoForcingLBGK<LatticeType>, LatticeType>() {}

        // Adds forcing to momentum
        void DoCalculateDensityMomentumFeq(
            HydroVars<GuoForcingLBGK>& hydroVars, site_t index);
        // Forwards to LBGK base class
        void DoCalculateFeq(HydroVars<GuoForcingLBGK>&, site_t);
        // Adds forcing to collision term
        void DoCollide(
          const LbmParameters* const lbmParams,
          HydroVars<GuoForcingLBGK>& hydroVars
        );
    };

  template<class LatticeType>
    struct HydroVars<GuoForcingLBGK<LatticeType> >
        : HydroVarsBase<LatticeType> {

        friend class GuoForcingLBGK<LatticeType>;
      public:
        // Pointer to force at this site
        const LatticeForceVector& force;

        template<class DataSource>
          HydroVars(geometry::Site<DataSource> const &_site)
            : HydroVarsBase<LatticeType>(_site), force(_site.GetForce()) {}

        HydroVars(
            const distribn_t* const f,
            const LatticeForceVector& _force
        ) : HydroVarsBase<LatticeType>(f), force(_force) {}

      protected:
        // Guo lattice distribution of external force contributions
        // as calculated in lattice::CalculateForceDistribution.
        FVector<LatticeType> forceDist;
    };

  template<class LatticeType>
    void GuoForcingLBGK<LatticeType> :: DoCalculateDensityMomentumFeq(
          HydroVars<GuoForcingLBGK<LatticeType> >& hydroVars, site_t index) {
      LatticeType::CalculateDensityMomentumFEq(hydroVars.f,
                                               hydroVars.force->x,
                                               hydroVars.force->y,
                                               hydroVars.force->z,
                                               hydroVars.density,
                                               hydroVars.momentum.x,
                                               hydroVars.momentum.y,
                                               hydroVars.momentum.z,
                                               hydroVars.velocity.x,
                                               hydroVars.velocity.y,
                                               hydroVars.velocity.z,
                                               hydroVars.f_eq.f);

      for (unsigned int ii = 0; ii < LatticeType::NUMVECTORS; ++ii)
        hydroVars.f_neq.f[ii] = hydroVars.f[ii] - hydroVars.f_eq.f[ii];
  }

  template<class LatticeType>
    void GuoForcingLBGK<LatticeType> :: DoCalculateFeq(
            HydroVars<GuoForcingLBGK>& hydroVars, site_t index) {
      LatticeType::CalculateFeq(hydroVars.density,
                                hydroVars.momentum.x,
                                hydroVars.momentum.y,
                                hydroVars.momentum.z,
                                hydroVars.f_eq.f);

      for (unsigned int ii = 0; ii < LatticeType::NUMVECTORS; ++ii)
        hydroVars.f_neq.f[ii] = hydroVars.f[ii] - hydroVars.f_eq.f[ii];
    }

  template<class LatticeType>
    void GuoForcingLBGK<LatticeType> :: DoCollide(
        const LbmParameters* const lbmParams,
        HydroVars<GuoForcingLBGK>& hydroVars
    ) {
      LatticeType::CalculateForceDistribution(
          lbmParams->GetTau(),
          hydroVars.velocity.x, hydroVars.velocity.y, hydroVars.velocity.z,
          hydroVars.force.x, hydroVars.force.y, hydroVars.force.z,
          hydroVars.forceDist.f
      );

      for (Direction dir(0); dir < LatticeType::NUMVECTORS; ++dir)
        hydroVars.SetFPostCollision(
            dir,
            hydroVars.f[dir]
            + hydroVars.f_neq.f[dir] * lbmParams->GetOmega()
            + hydroVars.forceDist.f[dir]
        );
  };

}}}

#endif /* HEMELB_LB_KERNELS_LBGK_H */
