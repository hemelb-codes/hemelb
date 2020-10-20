
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_TIMESTEP_TIMESTEPMANAGER_H
#define HEMELB_TIMESTEP_TIMESTEPMANAGER_H

#include <vector>

namespace hemelb
{
  namespace timestep
  {

    // A TimeStep is made up of one or more Phases of computation and
    // communication.
    
    // A Phase has multiple Actions which follow the same pattern:
    // 
    // * Begin - do any necessary computation/memory allocation before
    //   the receives begin
    // 
    // * Receive - start non-blocking receive operations
    //
    // * PreSend - do any necessary computation before data can be sent
    //
    // * Send - initiate non-blocking sends (or collectives)
    // 
    // * PreWait - do any available computation that does not depend
    //   on data to be received
    //
    // * Wait - wait for communication to complete
    //
    // * End - do computation that depends on received data

    // In addition there are two special Actions, BeginAll and EndAll
    // that will be executed before and after the main phases.

    // Thus for e.g. 2 phases, the ordering of Actions will be:
    // 
    // BeginAll0, BeginAll1
    // Begin0, Recv0, PreSend0, Send0, PreWait0, Wait0, End0,
    // Begin1, Recv1, PreSend1, Send1, PreWait1, Wait1, End1,
    // EndAll0, EndAll1

    // It is important to note that each Action can have multiple
    // Actors performing their tasks
    class Actor;
    
    // This class will arrange the calling of the different components
    // of HemeLB through a single timestep.
    class TimeStepManager
    {
    public:
      TimeStepManager(unsigned nPhases);
      // Register the actor for the phase
      //
      // NOTE: this is a non-owning pointer. You must manage the
      // lifetime of the actor. Needs a pointer to enable the
      // polymorphism.
      void AddToPhase(unsigned phase, Actor* actor);
      // Perform a timestep
      void DoStep();

    private:
      struct Phase
      {
	std::vector<Actor*> actors;
      };
      
      std::vector<Phase> mPhases;
    };
  }
}

#endif // HEMELB_TIMESTEP_TIMESTEPMANAGER_H
