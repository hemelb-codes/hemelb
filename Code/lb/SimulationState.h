// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

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
