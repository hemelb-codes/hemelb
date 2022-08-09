// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "net/ProcComms.h"
namespace hemelb
{
  namespace net
  {
    template <bool is_const>
    void ProcComms<is_const>::CreateMPIType()
    {
      std::vector<MPI_Aint> displacements(this->size());
      std::vector<int> lengths;
      std::vector<MPI_Datatype> types;

      int location = 0;

      MPI_Aint offset;
      MPI_Get_address(this->front().Pointer, &offset);

      for (auto& req: *this) {
        MPI_Get_address(req.Pointer, &displacements[location]);
        displacements[location] -= offset;

        ++location;
        lengths.push_back(req.Count);
        types.push_back(req.Type);
      }
      // Create the type and commit it.
      MPI_Type_create_struct(this->size(),
                             &lengths.front(),
                             &displacements.front(),
                             &types.front(),
                             &this->Type);
      MPI_Type_commit(&this->Type);
    }
  }
}
