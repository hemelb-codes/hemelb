// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "timestep/TimeStepManager.h"
#include "timestep/Actor.h"
#include "debug/Debugger.h"

namespace hemelb
{
  namespace timestep
  {
    TimeStepManager::TimeStepManager(unsigned nPhases) :
      mPhases(nPhases)
    {
    }
    // Register the actor for the phase
    void TimeStepManager::AddToPhase(unsigned phase, Actor* actor) {
      mPhases[phase].actors.push_back(actor);
    }
    
    void TimeStepManager::DoStep()
    {
      // BeginAll
      for (auto& ph: mPhases)
	for (auto ap: ph.actors)
	  ap->BeginAll();

      for (auto& ph: mPhases)
      {
	for (auto ap: ph.actors)
	  ap->Begin();
	
	for (auto ap: ph.actors)
	  ap->Receive();
	
	for (auto ap: ph.actors)
	  ap->PreSend();
	
	for (auto ap: ph.actors)
	  ap->Send();
	
	for (auto ap: ph.actors)
	  ap->PreWait();
	
	for (auto ap: ph.actors)
	  ap->Wait();
		
	for (auto ap: ph.actors)
	  ap->End();
      }

      for (auto& ph: mPhases)
	for (auto ap: ph.actors)
	  ap->EndAll();

    }

  }
}
