// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#include "lb/iolets/InOutLet.h"

namespace hemelb::lb::iolets
{
    void InOutLet::DoComms(const BoundaryCommunicator& bcComms, const LatticeTimeStep timeStep)
    {
        // pass
    }

    namespace {
        template<typename T>
        unsigned SmallestMagnitudeComponent(const util::Vector3D<T> &r) {
            // Return which direction has the smallest component
            const auto rSq = r.PointwiseMultiplication(r);
            return std::min_element(rSq.begin(), rSq.end()) - rSq.begin();
        }
    }

    // Dtor doesn't need to actually do anything.
    IoletExtraData::~IoletExtraData()
    {
    }

      IoletExtraData::IoletExtraData(InOutLet const& iolet) :
          n(iolet.GetNormal()), centre(iolet.GetPosition())
      {
        // Arbitrary vector - chosen as the direction mutually perpendicular
        // to the normal and a Cartesian direction that has the largest
        // magnitude.
        unsigned minInd = SmallestMagnitudeComponent(n);
        auto v = UnitVec::Zero();
        v[minInd] = 1.0;
        e1 = Cross(v, n).Normalise();
        e2 = Cross(n, e1);
      }

      LatticePosition IoletExtraData::WorldToIolet(LatticeVector r)
      {
        return WorldToIolet(LatticePosition{r});
      }

      LatticePosition IoletExtraData::WorldToIolet(LatticePosition r)
      {
        LatticePosition x = r - centre;
        return {Dot(e1, x), Dot(e2, x), Dot(n, x)};
      }

}
