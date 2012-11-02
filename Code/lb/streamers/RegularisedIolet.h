#ifndef HEMELB_LB_STREAMERS_REGULARISEDIOLET_H
#define HEMELB_LB_STREAMERS_REGULARISEDIOLET_H

#include "lb/streamers/BaseStreamer.h"
#include "util/utilityFunctions.h"
#include "debug/Debugger.h"

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
              geometry::Site<geometry::LatticeData> site = latDat->GetSite(siteIndex);

              const distribn_t* f = site.GetFOld<LatticeType> ();

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
              for (Direction direction = 1; direction < LatticeType::NUMVECTORS; ++direction)
              {
                if (!site.HasIolet(direction))
                  continue;

                int boundaryId = site.GetBoundaryId();

                // Set the density at the "ghost" site to be the density of the iolet.
                distribn_t ghostDensity = iolet->GetBoundaryDensity(boundaryId);

                // Calculate the velocity at the ghost site, as the component normal to the iolet.
                util::Vector3D<float> ioletNormal = iolet->GetLocalIolet(boundaryId)->GetNormal();

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

              BaseStreamer<RegularisedIolet>::template UpdateMinsAndMaxes<tDoRayTracing>(site,
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

#endif /* HEMELB_LB_STREAMERS_REGULARISEDIOLET_H */
