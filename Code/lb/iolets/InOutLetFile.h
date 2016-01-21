
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_IOLETS_INOUTLETFILE_H
#define HEMELB_LB_IOLETS_INOUTLETFILE_H

#include "lb/iolets/InOutLet.h"

namespace hemelb
{
  namespace lb
  {
    namespace iolets
    {

      /*
       * Template values are chosen to be tUpdatePeriod = 0 and tComms = false
       * If a trace is read from a file it should be done once and then stored
       * on each relevant proc. If memory is a concern tComms can be set to true
       * and then only the BCproc will keep the entire trace in memory
       * WARNING: - be cautious of setting tUpdatePeriod to something else other than
       * zero, because it may not be what you expect - see comments on CalculateCycle in
       * cc file.
       */
      class InOutLetFile : public InOutLet
      {
        public:
          InOutLetFile();
          virtual ~InOutLetFile();
          virtual InOutLet* Clone() const;
          virtual void Reset(SimulationState &state)
          {
            CalculateTable(state.GetTotalTimeSteps(), state.GetTimeStepLength());
          }

          const std::string& GetFilePath()
          {
            return pressureFilePath;
          }
          void SetFilePath(const std::string& path)
          {
            pressureFilePath = path;
          }

          LatticeDensity GetDensityMin() const
          {
            return densityMin;
          }
          LatticeDensity GetDensityMax() const
          {
            return densityMax;
          }
          LatticeDensity GetDensity(LatticeTimeStep timeStep) const
          {
            return densityTable[timeStep];
          }
          virtual void Initialise(const util::UnitConverter* unitConverter);
        private:
          void CalculateTable(LatticeTimeStep totalTimeSteps, PhysicalTime timeStepLength);
          std::vector<LatticeDensity> densityTable;
          LatticeDensity densityMin;
          LatticeDensity densityMax;
          std::string pressureFilePath;
          const util::UnitConverter* units;
      };

    }
  }
}

#endif /* HEMELB_LB_IOLETS_INOUTLETFILE_H */
