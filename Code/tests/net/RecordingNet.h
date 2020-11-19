// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_TESTS_NET_RECORDINGNET_H
#define HEMELB_TESTS_NET_RECORDINGNET_H
#include <map>
#include <cstdlib>
#include <iostream>

#include "constants.h"
#include "net/mpi.h"
#include "net/net.h"
#include "tests/net/LabelledRequest.h"

namespace hemelb
{
  namespace tests
  {
    namespace net
    {
      using namespace hemelb::net;
      /**
       * A mock for the net class for testing.
       * You must first give a complete, per-rank-ordered list of the sends and
       * another of the receives that it is to carry out.
       *
       * For sends, it checks that the sent data is identical to that specified/
       *
       * For recvs, it delivers the data that you specified.
       *
       * For both, it checks metadata.
       *
       * Once you have completed a "round" of communications, call
       * ExpectationsAllCompleted to check that there are none outstanding.
       *
       */
      class RecordingNet : public virtual StoringNet
      {
      public:
	RecordingNet(const MpiCommunicator& comms);

	/**
	 * Specify that this rank should receive a message
	 * @param pointer - to the message data that will later be received
	 * @param count - of number of elements in the message
	 * @param rank - of the source of the message
	 * @param label - used in error reporting (should be unique)
	 */
	template<class T>
	void RequireReceive(T* pointer, unsigned int count, proc_t rank,
			    const std::string &label = "")
	{
	  requiredReceipts[rank].push_back(LabelledRequest(pointer,
							   count,
							   MpiDataType<T>(),
							   rank,
							   label));
	}
	/**
	 * Specify that this rank should send a message
	 * @param pointer - to the message data that should be sent
	 * @param count - of the number of elements in the message
	 * @param rank - of the destination of the message
	 * @param label - used in error reporting (should be unique)
	 */
	template<class T>
	void RequireSend(T* pointer, unsigned int count, proc_t rank,
			 const std::string &label = "")
	{
	  requiredSends[rank].push_back(LabelledRequest(pointer,
							count,
							MpiDataType<T>(),
							rank,
							label));
	}

	/**
	 * Mock-execute queued receives.
	 *
	 * Does no actual communication, but copies in the mock data
	 * supplied to RequireReceive. It also checks that the receives
	 * match the required ones.
	 */
	void ReceivePointToPoint();
	
	/**
	 * Mock-execute queued sends.
	 *
	 * Does no actual communication, but checks that the sent data
	 * matches the mock data supplied to RequireSend. It also checks that
	 * the sends match the required ones.
	 */
	void SendPointToPoint();

	/**
	 * Mock-wait - clears the message queue
	 */
	void WaitPointToPoint();

	/**
	 * Assert that all required sends and receives have occurred.
	 */
	void ExpectationsAllCompleted();
      private:

	std::map<proc_t, BaseProcComms<LabelledRequest> > requiredReceipts;
	std::map<proc_t, BaseProcComms<LabelledRequest> > requiredSends;
      };

    }
  }
}
#endif
