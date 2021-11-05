// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_TESTS_NET_LABELLEDREQUEST_H
#define HEMELB_TESTS_NET_LABELLEDREQUEST_H

#include <string>

#include <mpi.h>

#include "net/StoredRequest.h"

namespace hemelb
{
  namespace tests
  {
    namespace net
    {
      class LabelledRequest : public hemelb::net::SimpleRequest {
      public:
	const std::string Label;
	LabelledRequest(void *pointer, int count, MPI_Datatype type, proc_t rank, const std::string &label);
	virtual bool EnvelopeIdentical(const SimpleRequest & other);
	virtual bool PayloadIdentical(const SimpleRequest & other);
	virtual void Unpack(SimpleRequest & other);
      };
    }
  }
}

#endif
