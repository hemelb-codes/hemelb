#ifndef HEMELB_LB_STREAMERS_NASHBB_H
#define HEMELB_LB_STREAMERS_NASHBB_H

#include "lb/streamers/BaseStreamer.h"
#include "util/utilityFunctions.h"

namespace hemelb
{
  namespace lb
  {
    namespace streamers
    {
      template<typename CollisionImpl>
      class NashBB : public BaseStreamer<NashBB<CollisionImpl> >
      {
        public:
          typedef CollisionImpl CollisionType;

        private:
          CollisionType collider;
          typedef typename CollisionType::CKernel::LatticeType LatticeType;

        public:
          NashBB(kernels::InitParams& initParams) :
            collider(initParams), iolet(initParams.boundaryObject)
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

              const distribn_t* f = site.GetFOld<LatticeType> ();

              kernels::HydroVars<typename CollisionType::CKernel> hydroVars(f);

              // First calculate the density and macro-velocity
              collider.CalculatePreCollision(hydroVars, site);

              collider.Collide(lbmParams, hydroVars);

              // Start by doing the normal stream and collide operation.
              for (Direction direction = 0; direction < LatticeType::NUMVECTORS; ++direction)
              {
                // The actual bounce-back lines, including streaming and collision. Basically swap
                // the non-equilibrium components of f in each of the opposing pairs of directions.
                site_t streamingDestination = site.HasBoundary(direction)
                  ? (siteIndex * LatticeType::NUMVECTORS) + LatticeType::INVERSEDIRECTIONS[direction]
                  : site.GetStreamedIndex<LatticeType> (direction);

                // Remember, oFNeq currently hold the equilibrium distribution. We
                // simultaneously use this and correct it, here.
                * (latDat->GetFNew(streamingDestination)) = hydroVars.GetFPostCollision()[direction];
              }

              // Let's next fill in the blanks on this site that won't get streamed to because of Iolets
              for (Direction direction = 1; direction < LatticeType::NUMVECTORS; ++direction)
              {
                if (!site.HasIolet(direction))
                  continue;

                distribn_t wall_distance = site.GetWallDistance<LatticeType> (direction);
                int boundaryId = site.GetBoundaryId();

                // Extrapolate from this site across the iolet to get the density at the
                // "ghost" site, which would be streaming to this one.
                distribn_t ghostDensity = (iolet->GetBoundaryDensity(boundaryId) + (wall_distance - 1.)
                    * hydroVars.density) / wall_distance;

                // Enforce a maximum 2% difference from the fluid site.
                ghostDensity = util::NumericalFunctions::enforceBounds(ghostDensity,
                                                                       hydroVars.density * 0.98,
                                                                       hydroVars.density * 1.02);

                // Calculate the velocity at the ghost site, as the component normal to the iolet.
                util::Vector3D<float> ioletNormal = iolet->GetLocalIolet(boundaryId)->GetNormal().GetNormalised();

                // Note that the division by density compensates for the fact that v_x etc have momentum
                // not velocity.
                distribn_t component = (hydroVars.momentum / hydroVars.density).Dot(ioletNormal);

                // TODO it's ugly that we have to do this.
                // TODO having to give 0 as an argument is also ugly.
                // TODO it's ugly that we have to give hydroVars a nonsense distribution vector
                // that doesn't get used.
                kernels::HydroVars<typename CollisionType::CKernel> ghostHydrovars(f);

                ghostHydrovars.density = ghostDensity;
                ghostHydrovars.momentum = ioletNormal * component * ghostDensity;

                collider.kernel.CalculateFeq(ghostHydrovars, 0);

                Direction unstreamed = LatticeType::INVERSEDIRECTIONS[direction];

                *latDat->GetFNew(siteIndex * LatticeType::NUMVECTORS + unstreamed)
                    = ghostHydrovars.GetFEq()[unstreamed];
              }

              ///< @todo #126 It would be nicer if tau is handled in a single place.
              hydroVars.tau = lbmParams->GetTau();

              BaseStreamer<NashBB>::template UpdateMinsAndMaxes<tDoRayTracing>(site,
                                                                               hydroVars,
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

        private:
          boundaries::BoundaryValues* iolet;
      };
    }
  }
}

#endif /* HEMELB_LB_STREAMERS_NASHBB_H */
