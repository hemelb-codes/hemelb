//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_LB_IOLETS_INOUTLETPARABOLICVELOCITY_H
#define HEMELB_LB_IOLETS_INOUTLETPARABOLICVELOCITY_H
#include "lb/iolets/InOutLetVelocity.h"

namespace hemelb
{
  namespace lb
  {
    namespace iolets
    {
      class InOutLetParabolicVelocity : public InOutLetVelocity
      {
        public:
          InOutLetParabolicVelocity();
          virtual ~InOutLetParabolicVelocity();
          void DoIO(TiXmlElement *iParent, bool iIsLoading, configuration::SimConfig* iSimConfig);
          InOutLet* Clone() const;
          LatticeVelocity GetVelocity(const LatticePosition& x, const LatticeTime t) const;
      };
    }
  }
}
#endif // HEMELB_LB_IOLETS_INOUTLETPARABOLICVELOCITY_H
