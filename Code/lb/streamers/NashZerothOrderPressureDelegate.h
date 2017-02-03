#ifndef HEMELB_LB_STREAMERS_NASHZEROTHORDERPRESSUREDELEGATE_H
#define HEMELB_LB_STREAMERS_NASHZEROTHORDERPRESSUREDELEGATE_H

#include "util/utilityFunctions.h"
#include "lb/streamers/BaseStreamerDelegate.h"

namespace hemelb
{
  namespace lb
  {
    namespace streamers
    {
      template<typename CollisionImpl>
      class NashZerothOrderPressureDelegate : public BaseStreamerDelegate<CollisionImpl>
      {
        public:
          typedef CollisionImpl CollisionType;
          typedef typename CollisionType::CKernel::LatticeType LatticeType;

          NashZerothOrderPressureDelegate(CollisionType& delegatorCollider,
                                          kernels::InitParams& initParams) :
              collider(delegatorCollider), iolet(*initParams.boundaryObject)
          {
          }

          inline void StreamLink(const LbmParameters* lbmParams,
                                 geometry::LatticeData* const latticeData,
                                 const geometry::Site<geometry::LatticeData>& site,
                                 kernels::HydroVars<typename CollisionType::CKernel>& hydroVars,
                                 const Direction& direction)
          {
            int boundaryId = site.GetIoletId();

            // Set the density at the "ghost" site to be the density of the iolet.
            distribn_t ghostDensity = iolet.GetBoundaryDensity(boundaryId);

            // Calculate the velocity at the ghost site, as the component normal to the iolet.
            util::Vector3D<float> ioletNormal = iolet.GetLocalIolet(boundaryId)->GetNormal();

            // Note that the division by density compensates for the fact that v_x etc have momentum
            // not velocity.
            distribn_t component = (hydroVars.momentum / hydroVars.density).Dot(ioletNormal);

            // TODO it's ugly that we have to do this.
            // TODO having to give 0 as an argument is also ugly.
            // TODO it's ugly that we have to give hydroVars a nonsense distribution vector
            // that doesn't get used.
            kernels::HydroVars<typename CollisionType::CKernel> ghostHydrovars(site);

            ghostHydrovars.density = ghostDensity;
            ghostHydrovars.momentum = ioletNormal * component * ghostDensity;

            collider.kernel.CalculateFeq(ghostHydrovars, 0);

            Direction unstreamed = LatticeType::INVERSEDIRECTIONS[direction];

            *latticeData->GetFNew(site.GetIndex() * LatticeType::NUMVECTORS + unstreamed) =
                ghostHydrovars.GetFEq()[unstreamed];
          }
        protected:
          CollisionType& collider;
          iolets::BoundaryValues& iolet;
      };
    }
  }
}

#endif // HEMELB_LB_STREAMERS_NASHZEROTHORDERPRESSUREDELEGATE_H
