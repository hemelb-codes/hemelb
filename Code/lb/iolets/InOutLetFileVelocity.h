//
// Copyright (C) University College London, 2007-2013, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_LB_IOLETS_INOUTLETFILEVELOCITY_H
#define HEMELB_LB_IOLETS_INOUTLETFILEVELOCITY_H

#include "InOutLetVelocity.h"

namespace hemelb
{
  namespace lb
  {
    namespace iolets
    {

      class InOutLetFileVelocity : public InOutLetVelocity
      {
        public:
          InOutLetFileVelocity();

          InOutLet* Clone() const;
          void Reset(SimulationState &state)
          {
            CalculateTable(state.GetTotalTimeSteps());
          }

          const std::string& GetFilePath()
          {
            return velocityFilePath;
          }
          void SetFilePath(const std::string& path)
          {
            velocityFilePath = path;
          }

          LatticeVelocity GetVelocity(const LatticePosition& x, const LatticeTimeStep t) const;

          void Initialise(const util::UnitConverter* unitConverter);

        private:
          std::string velocityFilePath;
          void CalculateTable(LatticeTimeStep totalTimeSteps);
          std::vector<LatticeSpeed> velocityTable;
          const util::UnitConverter* units;

      };

    }
  }
}
#endif /* HEMELB_LB_IOLETS_INOUTLETFILEVELOCITY_H */
