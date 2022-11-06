// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_TESTS_HELPERS_LABELLEDREQUEST_H
#define HEMELB_TESTS_HELPERS_LABELLEDREQUEST_H

#include <string>

#include <mpi.h>

#include "net/StoredRequest.h"

namespace hemelb
{
  namespace tests
  {
    namespace net
    {
        template <bool is_const>
    class LabelledRequest : public hemelb::net::SimpleRequest<is_const> {
      public:
        using Base = hemelb::net::SimpleRequest<is_const>;
	const std::string Label;
	LabelledRequest(typename Base::ptr pointer, int count, MPI_Datatype type, proc_t rank, const std::string &label);
	bool EnvelopeIdentical(const Base& other);
	bool PayloadIdentical(const Base & other);
	void Unpack(Base& other);
      };

        // Declare that we will explicitly instantiate
        extern template class LabelledRequest<true>;
        extern template class LabelledRequest<false>;
    }
  }
}

#endif
