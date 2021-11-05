// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <catch2/catch.hpp>

#include "tests/net/RecordingNet.h"

namespace hemelb
{
  namespace tests
  {
    namespace net
    {
      RecordingNet::RecordingNet(const MpiCommunicator& comms) :
	BaseNet(comms), StoringNet(comms), requiredReceipts(), requiredSends()
      {
      }

      /**
       * Mock-execute queued receives.
       *
       * Does no actual communication, but copies in the mock data
       * supplied to RequireReceive. It also checks that the receives
       * match the required ones.
       */
      void RecordingNet::ReceivePointToPoint()
      {
	for (std::map<proc_t, ProcComms>::iterator it = receiveProcessorComms.begin();
	     it != receiveProcessorComms.end(); ++it)
	  {
	    for (ProcComms::iterator message = it->second.begin(); message != it->second.end();
		 message++)
              {

                if (requiredReceipts[message->Rank].size() == 0)
		  {
		    // It should not be the case, that a request is made, but there is no recorded assertion in the list of recorded assertions.
		    // i.e., we've popped off all the recorded requests, and there are none left to match the request we just received.
		    REQUIRE(requiredReceipts[message->Rank].size() != 0);
		    return;
		  }

                REQUIRE(requiredReceipts[message->Rank].front().EnvelopeIdentical(*message));
                requiredReceipts[message->Rank].front().Unpack(*message);

                // This assertion just checks that "Unpack" has done it's job.
                // It shouldn't fail due to problems in tested code
                // Would only fail due to memory screwups or bad code in this class's Unpack method.
                REQUIRE(requiredReceipts[message->Rank].front().PayloadIdentical(*message));

                requiredReceipts[message->Rank].pop_front();
              }

	  }

      }
      /**
       * Mock-execute queued sends.
       *
       * Does no actual communication, but checks that the sent data
       * matches the mock data supplied to RequireSend. It also checks that
       * the sends match the required ones.
       */
      void RecordingNet::SendPointToPoint()
      {

	for (std::map<proc_t, ProcComms>::iterator it = sendProcessorComms.begin();
	     it != sendProcessorComms.end(); ++it)
	  {
	    for (ProcComms::iterator message = it->second.begin(); message != it->second.end();
		 message++)
              {

                if (requiredSends[message->Rank].size() == 0)
		  {
		    REQUIRE(requiredSends[message->Rank].size() != 0);
		    return;
		  }
                REQUIRE(requiredSends[message->Rank].front().EnvelopeIdentical(*message));
                REQUIRE(requiredSends[message->Rank].front().PayloadIdentical(*message));
                requiredSends[message->Rank].pop_front();
              }

	  }
      }

      /**
       * Mock-wait - clears the message queue
       */
      void RecordingNet::WaitPointToPoint()
      {
	receiveProcessorComms.clear();
	sendProcessorComms.clear();
      }

      /**
       * Assert that all required sends and receives have occurred.
       */
      void RecordingNet::ExpectationsAllCompleted()
      {
	for (std::map<proc_t, BaseProcComms<LabelledRequest> >::iterator receipts_from_core =
	       requiredReceipts.begin(); receipts_from_core != requiredReceipts.end();
	     receipts_from_core++)
	  {
	    REQUIRE(0 == receipts_from_core->second.size());
	  }
	for (std::map<proc_t, BaseProcComms<LabelledRequest> >::iterator sends_to_core =
	       requiredSends.begin(); sends_to_core != requiredSends.end(); sends_to_core++)
	  {
	    REQUIRE(0 == sends_to_core->second.size());
	  }
      }

    }
  }
}
