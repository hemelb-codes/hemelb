// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

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
            CalculateTable(state.GetTotalTimeSteps());
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
        private:
          void CalculateTable(LatticeTimeStep totalTimeSteps);
          std::vector<LatticeDensity> densityTable;
          LatticeDensity densityMin;
          LatticeDensity densityMax;
          std::string pressureFilePath;
      };

    }
  }
}

#endif /* HEMELB_LB_IOLETS_INOUTLETFILE_H */
