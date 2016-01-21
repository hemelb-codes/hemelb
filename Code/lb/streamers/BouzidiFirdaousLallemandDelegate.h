
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_STREAMERS_BOUZIDIFIRDAOUSLALLEMANDDELEGATE_H
#define HEMELB_LB_STREAMERS_BOUZIDIFIRDAOUSLALLEMANDDELEGATE_H

#include "lb/streamers/BaseStreamerDelegate.h"
#include "lb/streamers/SimpleBounceBackDelegate.h"

namespace hemelb
{
  namespace lb
  {
    namespace streamers
    {
      /**
       * Implement the boundary condition described in Bouzidi, Firdaous and
       * Lallemand "Momentum transfer of a Boltzmann-lattice fluid with
       * boundaries", Phys. Fluids 13, 3452-3459 (2001).
       *
       * This is based on the idea of doing interpolated bounce-back.
       *
       * Note that since the method requires data from neighbouring sites (in
       * some circumstances), it has a DoPostStep method.
       */
      template<typename CollisionImpl>
      class BouzidiFirdaousLallemandDelegate : public BaseStreamerDelegate<CollisionImpl>
      {
        public:
          typedef CollisionImpl CollisionType;
          typedef typename CollisionType::CKernel::LatticeType LatticeType;
        private:
          SimpleBounceBackDelegate<CollisionType> bbDelegate;

        public:
          BouzidiFirdaousLallemandDelegate(CollisionType& delegatorCollider, kernels::InitParams& initParams) :
            bbDelegate(delegatorCollider, initParams)
          {
          }

          inline void StreamLink(const LbmParameters* lbmParams,
                                 geometry::LatticeData* const latticeData,
                                 const geometry::Site<geometry::LatticeData>& site,
                                 kernels::HydroVars<typename CollisionType::CKernel>& hydroVars,
                                 const Direction& direction)
          {
            site_t invDirection = LatticeType::INVERSEDIRECTIONS[direction];
            site_t bbDestination = (site.GetIndex() * LatticeType::NUMVECTORS) + invDirection;
            distribn_t q = site.GetWallDistance<LatticeType> (direction);

            if (site.HasWall(invDirection) || q < 0.5)
            {
              // If there IS NO fluid site in the opposite direction, fall back to SBB.
              // If there IS such a site, we have to wait for the site in the opposite
              // direction to finish in order to complete this update. So just bounce-back
              // the post collision f that we would have otherwise thrown away (to avoid
              // having to collide twice).
              bbDelegate.StreamLink(lbmParams, latticeData, site, hydroVars, direction);
            }
            else
            {
              // We have a fluid site and have all the data needed to complete this direction!
              // Implement Eq (5b) from Bouzidi et al.
              * (latticeData->GetFNew(bbDestination)) = (hydroVars.GetFPostCollision()[direction] + (2.0 * q - 1)
                  * hydroVars.GetFPostCollision()[invDirection]) / (2.0 * q);
            }

          }
          inline void PostStepLink(geometry::LatticeData* const latticeData,
                                   const geometry::Site<geometry::LatticeData>& site,
                                   const Direction& direction)
          {
            distribn_t* fNew = latticeData->GetFNew(site.GetIndex() * LatticeType::NUMVECTORS);
            site_t invDirection = LatticeType::INVERSEDIRECTIONS[direction];
            distribn_t q = site.GetWallDistance<LatticeType> (direction);
            // If there is no fluid site in the opposite direction, fall back to simple
            // bounce back, which has been done above.

            // If q >= 0.5, then we handled that fully above also.
            if (!site.HasWall(invDirection) && q < 0.5)
            {
              // So, we have a fluid site and all the data needed to complete this direction!
              // Implement Eq (5a) from Bouzidi et al.

              // Note that:
              // - fNew[direction] is the newly-arrived fPostColl[direction] from the neighbouring site
              // - fNew[invDirection] is the above-bounced-back fPostColl[direction] for this site.
              fNew[invDirection] = 2.0 * q * fNew[invDirection] + (1.0 - 2.0 * q) * fNew[direction];
            }
          }
      };
    }
  }
}

#endif /* HEMELB_LB_STREAMERS_BOUZIDIFIRDAOUSLALLEMANDDELEGATE_H */
