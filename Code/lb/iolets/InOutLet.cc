
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#include "lb/iolets/InOutLet.h"

namespace hemelb
{
  namespace lb
  {
    namespace iolets
    {
      void InOutLet::DoComms(const BoundaryCommunicator& bcComms, const LatticeTimeStep timeStep)
      {
        // pass
      }

      namespace
      {
        unsigned SmallestMagnitudeComponent(const LatticeVector r)
        {
          // Return which direction has the smallest component
          LatticeVector rSq = r.PointwiseMultiplication(r);
          if (rSq.x < rSq.y)
          {
            if (rSq.x < rSq.z)
            {
              return 0;
            }
            else
            {
              return 2;
            }
          }
          else
          {
            if (rSq.y < rSq.z)
            {
              return 1;
            }
            else
            {
              return 2;
            }
          }
        }

      }

      // Dtor doesn't need to actually do anything.
      IoletExtraData::~IoletExtraData()
      {
      }

      IoletExtraData::IoletExtraData(InOutLet& iolet) :
        n(iolet.GetNormal()),
            centre(iolet.GetPosition())
      {
        // Arbitrary vector - chosen as the direction mutually perpendicular
        // to the normal and a Cartesian direction that has the largest
        // magnitude.
        unsigned minInd = SmallestMagnitudeComponent(n);
        util::Vector3D < Dimensionless > v(0);
        v[minInd] = 1.;
        e1 = UnitVec::Cross(v, n).Normalise();
        e2 = UnitVec::Cross(n, e1);
      }

      LatticePosition IoletExtraData::WorldToIolet(LatticeVector r)
      {
        LatticePosition x(r);
        x -= centre;
        return LatticePosition(e1.Dot(x), e2.Dot(x), n.Dot(x));
      }

      LatticePosition IoletExtraData::WorldToIolet(LatticePosition r)
      {
        LatticePosition x = r - centre;
        return LatticePosition(e1.Dot(x), e2.Dot(x), n.Dot(x));
      }

    }
  }
}
