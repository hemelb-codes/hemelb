
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#ifndef HEMELB_LB_STENTS_STENTFILEFLUX_H
#define HEMELB_LB_STENTS_STENTFILEFLUX_H

#include "lb/stents/StentFlux.h"

namespace hemelb
{
  namespace lb
  {
    namespace stents
    {

      class StentFileFlux : public StentFlux
      {
        public:
          StentFileFlux();

          Stent* Clone() const;
          void Reset(SimulationState &state)
          {
            CalculateTable(state.GetTotalTimeSteps(), state.GetTimeStepLength());
          }

          const std::string& GetFilePath()
          {
            return fluxFilePath;
          }
          void SetFilePath(const std::string& path)
          {
            fluxFilePath = path;
          }

          LatticeSpeed GetFlux(const LatticeTimeStep t) const;
          /*LatticeVelocity GetVelocity2(const util::Vector3D<int64_t> globalCoordinates,
                                                                  const LatticeTimeStep t) const;*/

          void Initialise(const util::UnitConverter* unitConverter);

          bool useWeightsFromFile;

        private:
          std::string fluxFilePath;
          std::string fluxWeightsFilePath;
          void CalculateTable(LatticeTimeStep totalTimeSteps, PhysicalTime timeStepLength);
          std::vector<LatticeSpeed> fluxTable;
          const util::UnitConverter* units;

          std::map<std::vector<int>, double> weights_table;

          //double calcVTot(std::vector<double> v);

          //std::vector<double> updateV(std::vector<double> v, std::vector<int> xyz, std::map<std::vector<int>, double> weights_table);

      };

    }
  }
}
#endif /* HEMELB_LB_STENTS_STENTFILEFLUX_H */
