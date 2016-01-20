
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_SIMULATIONSTATE_H
#define HEMELB_LB_SIMULATIONSTATE_H

#include "units.h"
#include "reporting/Reportable.h"
#include "units.h"

namespace hemelb
{
  namespace lb
  {
    enum Stability
    {
      UndefinedStability = -1,
      Unstable = 0,
      Stable = 1,
      StableAndConverged = 2
    };

    class SimulationState : public reporting::Reportable
    {
      public:
        SimulationState(double timeStepLength, unsigned long totalTimeSteps);

        void Increment();
        void Reset();
        void SetIsTerminating(bool value);
        void SetIsRendering(bool value);
        void SetStability(Stability value);

        LatticeTimeStep GetTimeStep() const;
        LatticeTimeStep Get0IndexedTimeStep() const;
        LatticeTimeStep GetTotalTimeSteps() const;
        bool IsTerminating() const;
        bool IsRendering() const;
        Stability GetStability() const;

        PhysicalTime GetTime() const {return GetTimeStepLength()*Get0IndexedTimeStep();}
        PhysicalTime GetTimeStepLength() const {return timeStepLength;}

        void Report(ctemplate::TemplateDictionary& dictionary);

      private:
        PhysicalTime timeStepLength;
        LatticeTimeStep timeStep;
        LatticeTimeStep totalTimeSteps;
        bool isTerminating;
        bool isRendering;
        Stability stability;
    };
  }
}

#endif /* SIMULATIONSTATE_H_ */
