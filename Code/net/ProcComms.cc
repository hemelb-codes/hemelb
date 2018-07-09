
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "net/ProcComms.h"
namespace hemelb
{
  namespace net
  {
    void ProcComms::CreateMPIType()
    {
      std::vector<MPI_Aint> displacements(size());
      std::vector<int> lengths;
      std::vector<MPI_Datatype> types;

      int location = 0;

      MPI_Aint offset;
      MPI_Get_address(front().Pointer, &offset);

      for (iterator it = begin(); it != end(); ++it)
      {
        MPI_Get_address(it->Pointer, &displacements[location]);
        displacements[location] -= offset;

        ++location;
        lengths.push_back(it->Count);
        types.push_back(it->Type);
      }
      // Create the type and commit it.
      MPI_Type_create_struct(this->size(), &lengths.front(), &displacements.front(), &types.front(), &Type);
      MPI_Type_commit(&Type);
    }
  }
}
