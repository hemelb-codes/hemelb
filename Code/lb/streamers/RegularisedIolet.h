#ifndef HEMELB_LB_STREAMERS_REGULARISEDIOLET_H
#define HEMELB_LB_STREAMERS_REGULARISEDIOLET_H

#include "lb/streamers/BaseStreamer.h"

namespace hemelb
{
  namespace lb
  {
    namespace streamers
    {
      template<typename CollisionImpl>
      class RegularisedIolet : public BaseStreamer<RegularisedIolet<CollisionImpl> >
      {
        public:
          typedef CollisionImpl CollisionType;

        private:
          CollisionType collider;
          typedef typename CollisionType::CKernel::LatticeType LatticeType;

        public:
          RegularisedIolet(kernels::InitParams& initParams) :
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

              distribn_t* f = site.GetFOld<LatticeType> ();

              kernels::HydroVars<typename CollisionType::CKernel> hydroVars(f);

              // First calculate the density and macro-velocity
              collider.CalculatePreCollision(hydroVars, site);

              collider.Collide(lbmParams, hydroVars);

              // Start by doing the normal stream and collide operation.
              for (Direction direction = 0; direction < LatticeType::NUMVECTORS; ++direction)
              {
                * (latDat->GetFNew(site.GetStreamedIndex<LatticeType> (direction)))
                    = hydroVars.GetFPostCollision()[direction];
              }

              // Let's next fill in the blanks on this site that won't get streamed to.
              for (Direction direction = 0; direction < LatticeType::NUMVECTORS; ++direction)
              {
                if (!site.HasIolet(direction))
                  continue;

                distribn_t wall_distance = site.GetWallDistance<LatticeType> (direction);
                int boundaryId = site.GetBoundaryId();

                // Extrapolate from this site across the iolet to get the density at the
                // "ghost" site, which would be streaming to this one.
                distribn_t ghostDensity = (iolet->GetBoundaryDensity(boundaryId) + (wall_distance - 1.)
                    * hydroVars.density) / wall_distance;

                // Calculate the velocity at the ghost site, as the component normal to the iolet.
                util::Vector3D<float> ioletNormal = iolet->GetLocalIolet(boundaryId)->GetNormal().GetNormalised();

                // Note that the division by density compensates for the fact that v_x etc have momentum
                // not velocity.
                float component = ioletNormal.Dot(util::Vector3D<float>(hydroVars.v_x, hydroVars.v_y, hydroVars.v_z))
                    / hydroVars.density;

                util::Vector3D<distribn_t> ghostVelocity = ioletNormal * component * ghostDensity;

                // TODO it's ugly that we have to do this.
                // TODO having to give 0 as an argument is also ugly.
                // TODO it's ugly that we have to give hydroVars a nonsense distribution vector
                // that doesn't get used.
                kernels::HydroVars<typename CollisionType::CKernel> ghostHydrovars(f);

                ghostHydrovars.density = ghostDensity;
                ghostHydrovars.v_x = ghostVelocity.x;
                ghostHydrovars.v_y = ghostVelocity.y;
                ghostHydrovars.v_z = ghostVelocity.z;

                collider.kernel.CalculateFeq(ghostHydrovars, 0);

                *latDat->GetFNew(siteIndex * LatticeType::NUMVECTORS + direction) = ghostHydrovars.GetFEq()[direction];
              }

              ///< @todo #126 It would be nicer if tau is handled in a single place.
              hydroVars.tau = lbmParams->GetTau();

              BaseStreamer<RegularisedIolet>::template UpdateMinsAndMaxes<tDoRayTracing>(hydroVars.v_x,
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

        private:
          boundaries::BoundaryValues* iolet;
      };
    }
  }
}

#endif /* HEMELB_LB_STREAMERS_REGULARISED_H */
