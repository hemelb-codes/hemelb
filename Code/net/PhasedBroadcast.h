#ifndef HEMELB_NET_PHASEDBROADCAST_H
#define HEMELB_NET_PHASEDBROADCAST_H

#include <vector>
#include "net/net.h"
#include "lb/SimulationState.h"
#include "topology/NetworkTopology.h"

namespace hemelb
{
  namespace net
  {
    class PhasedBroadcast
    {
      public:
        PhasedBroadcast(Net * iNet,
                        const lb::SimulationState * iSimState,
                        unsigned int spreadFactor);

        virtual void Reset();

        void RequestComms();
        void PostReceive();

      protected:
        virtual void ProgressFromChildren();
        virtual void ProgressFromParent();
        virtual void ProgressToChildren();
        virtual void ProgressToParent();

        virtual void Effect();

        template<class T>
        void SendToChildren(T* data, int count);

        template<class T>
        void ReceiveFromParent(T* data, int count);

        /**
         * Receives data from each child. This is a set length per child, and each child's
         * data is inserted contiguously into the provided array.
         *
         * @param dataStart Pointer to the start of the array.
         * @param countPerChild Number of elements to receive per child.
         */
        template<class T>
        void ReceiveFromChildren(T* dataStart, int countPerChild);

        template<class T>
        void SendToParent(T* data, int count);

      private:
        static const int NOPARENT = -1;

        int Get0IndexedCycleNumber() const;

        /**
         * Note that depths are 0-indexed. I.e. a tree with a single node has
         * a depth of zero, and the node itself is considered to be at depth 0.
         */
        unsigned int mMyDepth;
        unsigned int mTreeDepth;
        unsigned int mSpreadFactor;
        int mParent;
        std::vector<int> mChildren;

        const lb::SimulationState * mSimState;
        Net * mNet;
        topology::NetworkTopology * mNetTop;
    };
  }
}

#endif /* HEMELB_NET_PHASEDBROADCAST_H */
