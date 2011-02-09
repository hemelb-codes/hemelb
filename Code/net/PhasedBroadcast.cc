#include <math.h>

#include "net/PhasedBroadcast.h"
#include "net/net.h"
#include "lb/LbmParameters.h"

namespace hemelb
{
  namespace net
  {
    PhasedBroadcast::PhasedBroadcast(Net* iNet,
                                     const lb::SimulationState * iSimState,
                                     unsigned int spreadFactor)
    {
      mNet = iNet;
      mSimState = iSimState;

      mTreeDepth = 0;
      mMyDepth = 0;
      mSpreadFactor = spreadFactor;
      int noSeenToThisDepth = 1;
      int noAtCurrentDepth = 1;

      while (noSeenToThisDepth < mNetTop->GetProcessorCount())
      {
        ++mTreeDepth;
        noAtCurrentDepth *= spreadFactor;
        noSeenToThisDepth += noAtCurrentDepth;

        if (noSeenToThisDepth >= mNetTop->GetLocalRank() && (noSeenToThisDepth - noAtCurrentDepth)
            < mNetTop->GetLocalRank())
        {
          mMyDepth = mTreeDepth;
        }
      }

      // In a M-tree, with a root of 0, each node N's parent has rank floor((N-1) / M)
      if (mNetTop->GetLocalRank() == 0)
      {
        mParent = NOPARENT;
      }
      else
      {
        mParent = (mNetTop->GetLocalRank() - 1) / spreadFactor;
      }

      // The children of a node N in a M-tree with root 0 are those in the range (M*N)+1,...,(M*N) + M
      for (unsigned int child = (spreadFactor * mNetTop->GetLocalRank()) + 1; child
          <= (spreadFactor + 1) * mNetTop->GetLocalRank(); ++child)
      {
        mChildren.push_back(child);
      }
    }

    /*
     * PROGRAMME
     * ---------
     *
     * For 0-indexed cycle counts and 0-indexed depths, tree depth = N.
     *
     * Iterations 0 to (N-1) are nodes at depth (it) passing down to nodes at depth (it + 1).
     * Action happens on all nodes at the end of iteration (N-1).
     * Iterations N to (2N-1) are nodes at depth (2N - it) passing up to nodes at depth (2N - it - 1).
     */
    void PhasedBroadcast::RequestComms()
    {
      int iCycleNumber = Get0IndexedCycleNumber();

      if (mMyDepth == iCycleNumber)
      {
        ProgressToChildren();
      }
      else if (mMyDepth == (iCycleNumber + 1))
      {
        ProgressFromParent();
      }
      else if (mMyDepth == (2 * mTreeDepth - iCycleNumber))
      {
        ProgressToParent();
      }
      else if (mMyDepth == (2 * mTreeDepth - iCycleNumber - 1))
      {
        ProgressFromChildren();
      }
    }

    void PhasedBroadcast::PostReceive()
    {
      if (Get0IndexedCycleNumber() == (mTreeDepth - 1))
      {
        Effect();
      }
    }

    int PhasedBroadcast::Get0IndexedCycleNumber() const
    {
      return (mSimState->CycleId - 1) % (2 * mTreeDepth);
    }

    /**
     * The default version of this method does nothing.
     */
    void PhasedBroadcast::Reset()
    {

    }

    /**
     * The default version of this method does nothing.
     */
    void PhasedBroadcast::ProgressFromChildren()
    {

    }

    /**
     * The default version of this method does nothing.
     */
    void PhasedBroadcast::ProgressFromParent()
    {

    }

    /**
     * The default version of this method does nothing.
     */
    void PhasedBroadcast::ProgressToChildren()
    {

    }

    /**
     * The default version of this method does nothing.
     */
    void PhasedBroadcast::ProgressToParent()
    {

    }

    /**
     * The default version of this method does nothing.
     */
    void PhasedBroadcast::Effect()
    {

    }

    template<class T>
    void PhasedBroadcast::SendToChildren(T* data, int count)
    {
      for (std::vector<int>::iterator it = mChildren.begin(); it != mChildren.end(); ++it)
      {
        mNet->RequestSend<T> (data, count, *it);
      }
    }

    template<class T>
    void PhasedBroadcast::ReceiveFromParent(T* data, int count)
    {
      if (mParent != NOPARENT)
      {
        mNet ->RequestReceive<T> (data, count, mParent);
      }
    }

    /**
     * Receives data from each child. This is a set length per child, and each child's
     * data is inserted contiguously into the provided array.
     *
     * @param dataStart Pointer to the start of the array.
     * @param countPerChild Number of elements to receive per child.
     */
    template<class T>
    void PhasedBroadcast::ReceiveFromChildren(T* dataStart, int countPerChild)
    {
      T* data = dataStart;
      for (std::vector<int>::iterator it = mChildren.begin(); it != mChildren.end(); ++it)
      {
        mNet->RequestReceive<T> (data, countPerChild, *it);
        data += countPerChild;
      }
    }

    template<class T>
    void PhasedBroadcast::SendToParent(T* data, int count)
    {
      if (mParent != NOPARENT)
      {
        mNet->RequestSend<T> (data, count, mParent);
      }
    }
  }
}
