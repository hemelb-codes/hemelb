#ifndef HEMELB_NET_ITERATEDACTION_H
#define HEMELB_NET_ITERATEDACTION_H

namespace hemelb
{
  namespace net
  {
    class IteratedAction
    {
      public:
        virtual void RequestComms();
        virtual void PreSend();
        virtual void PreReceive();
        virtual void PostReceive();
        virtual void EndIteration();
        virtual void Reset();
    };
  }
}

#endif /* HEMELB_NET_ITERATEDACTION_H */
