// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_NET_ITERATEDACTION_H
#define HEMELB_NET_ITERATEDACTION_H
#include "net/phased/Concern.h"
namespace hemelb
{
  namespace net
  {
    class IteratedAction : public phased::Concern
    {
      public:
      /***
       * When an iterated actor is called through the phased::StepManager mechanism,
       * this method dispatches to the individual step methods
       * @param action Enumeration indicating the step
       * @return True if an action was successfully called for the step
       */
        bool CallAction(int action);
        virtual ~IteratedAction();
        virtual void RequestComms();
        virtual void PreSend();
        virtual void PreReceive();
        virtual void PostReceive();
        virtual void EndIteration();
     };
  }
}

#endif /* HEMELB_NET_ITERATEDACTION_H */
