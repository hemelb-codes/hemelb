//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#include "comm/MpiRequest.h"

namespace hemelb
{
  namespace comm
  {
    MpiRequest::MpiRequest() :
        req(MPI_REQUEST_NULL)
    {
    }
    MpiRequest::MpiRequest(MPI_Request req_) :
        req(req_)
    {
    }

    MpiRequest::operator bool() const
    {
      return req != MPI_REQUEST_NULL;
    }

    void MpiRequest::Wait()
    {
      HEMELB_MPI_CALL(MPI_Wait, (&req, MPI_STATUS_IGNORE));
    }

    bool MpiRequest::Test()
    {
      int flag;
      HEMELB_MPI_CALL(MPI_Test, (&req, &flag, MPI_STATUS_IGNORE));
      return flag;
    }
    
    size_t MpiRequestList::size() const
    {
      return reqs.size();
    }
    
    void MpiRequestList::resize(size_t i)
    {
      reqs.resize(i);
    }
    
    void MpiRequestList::push_back(Request::Ptr r)
    {
      auto mr = std::dynamic_pointer_cast<MpiRequest>(r);
      //MpiRequest&& mr = dynamic_cast<MpiRequest&&>(r);
      reqs.push_back(mr->req);
      mr->req = MPI_REQUEST_NULL;
    }
    
    void MpiRequestList::set(size_t i, Request::Ptr r) {
      auto mr = std::dynamic_pointer_cast<MpiRequest>(r);
      //MpiRequest&& mr = dynamic_cast<MpiRequest&&>(r);
      reqs[i] = mr->req;
      mr->req = MPI_REQUEST_NULL;      
    }
    
    void MpiRequestList::WaitAll()
    {
      HEMELB_MPI_CALL(
		      MPI_Waitall,
		      (reqs.size(), reqs.data(), MPI_STATUSES_IGNORE)
		      );
    }
    
    bool MpiRequestList::TestAll()
    {
      int flag;
      HEMELB_MPI_CALL(
		      MPI_Testall,
		      (reqs.size(), reqs.data(), &flag, MPI_STATUSES_IGNORE)
		      );
      return flag;
    }
    
  }
}
