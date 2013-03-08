//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_LB_BOUNDARIES_IOLETS_INOUTLETPARABOLICVELOCITY_H
#define HEMELB_LB_BOUNDARIES_IOLETS_INOUTLETPARABOLICVELOCITY_H
#include "lb/boundaries/iolets/InOutLet.h"

namespace hemelb
{
  namespace lb
  {
    namespace boundaries
    {
      namespace iolets
      {
        class InOutLetParabolicVelocity : public InOutLet
        {
          public:
            InOutLetParabolicVelocity();
            virtual ~InOutLetParabolicVelocity();
            virtual void DoIO(TiXmlElement *iParent, bool iIsLoading, configuration::SimConfig* iSimConfig);
            virtual InOutLet* Clone() const;
            virtual PhysicalPressure GetPressureMin() const;
            virtual PhysicalPressure GetPressureMax() const;
            virtual LatticeDensity GetDensity(LatticeTime time_step) const;

            virtual void Reset(SimulationState &state)
            {
              //pass;
            }
            /**
             * Note that the radius and max speed for these are specified in LATTICE UNITS in the XML file.
             * This is indeed a horrible hack.
             * @return
             */
            virtual LatticeDistance& GetRadius()
            {
              return radius;
            }
            virtual void SetRadius(LatticeDistance r)
            {
              radius = r;
            }

            virtual LatticeSpeed& GetMaxSpeed()
            {
              return maxSpeed;
            }
            virtual void SetMaxSpeed(LatticeSpeed v)
            {
              maxSpeed = v;
            }


            virtual LatticeVelocity GetVelocityAtPosition(LatticePosition x);

          protected:
            LatticeDistance radius;
            LatticeSpeed maxSpeed;
        };
      }
    }
  }
}
#endif // HEMELB_LB_BOUNDARIES_IOLETS_INOUTLETPARABOLICVELOCITY_H
