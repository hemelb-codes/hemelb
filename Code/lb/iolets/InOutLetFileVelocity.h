
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#ifndef HEMELB_LB_IOLETS_INOUTLETFILEVELOCITY_H
#define HEMELB_LB_IOLETS_INOUTLETFILEVELOCITY_H

#include "lb/iolets/InOutLetVelocity.h"

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
            CalculateTable(state.GetTotalTimeSteps(), state.GetTimeStepLength());
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
          /*LatticeVelocity GetVelocity2(const util::Vector3D<int64_t> globalCoordinates,
                                                                  const LatticeTimeStep t) const;*/

          void Initialise(const util::UnitConverter* unitConverter);

          bool useWeightsFromFile;

        private:
          std::string velocityFilePath;
          std::string velocityWeightsFilePath;
          void CalculateTable(LatticeTimeStep totalTimeSteps, PhysicalTime timeStepLength);
          std::vector<LatticeSpeed> velocityTable;
          const util::UnitConverter* units;

          std::map<std::vector<int>, double> weights_table;

          //double calcVTot(std::vector<double> v);

          //std::vector<double> updateV(std::vector<double> v, std::vector<int> xyz, std::map<std::vector<int>, double> weights_table);

      };

    }
  }
}
#endif /* HEMELB_LB_IOLETS_INOUTLETFILEVELOCITY_H */
