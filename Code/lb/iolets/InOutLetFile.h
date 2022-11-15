// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_IOLETS_INOUTLETFILE_H
#define HEMELB_LB_IOLETS_INOUTLETFILE_H

#include <filesystem>
#include <utility>

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
          virtual ~InOutLetFile() override = default;
          InOutLet* clone() const override;
          void Reset(SimulationState& state) override;

          inline const std::filesystem::path& GetFilePath()
          {
            return pressureFilePath;
          }
          inline void SetFilePath(const std::string& path)
          {
            pressureFilePath = path;
          }

          inline LatticeDensity GetDensityMin() const override
          {
            return densityMin;
          }
          inline LatticeDensity GetDensityMax() const override
          {
            return densityMax;
          }
          inline LatticeDensity GetDensity(LatticeTimeStep timeStep) const override
          {
            return densityTable[timeStep];
          }
          void Initialise(const util::UnitConverter* unitConverter) override;

        private:
          std::vector<LatticeDensity> densityTable;
          LatticeDensity densityMin;
          LatticeDensity densityMax;
          std::filesystem::path pressureFilePath;
          using DataPair = std::pair<LatticeTime, LatticeDensity>;
          std::vector<DataPair> file_data_lat;
      };

    }
  }
}

#endif /* HEMELB_LB_IOLETS_INOUTLETFILE_H */
