//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#include "geometry/GeometryPeeker.h"
#include "net/NetworkTopology.h"

namespace hemelb
{
  namespace geometry
  {

    GeometryPeeker::GeometryPeeker(const std::string& gmyFileName)
    {

      currentComms = net::NetworkTopology::Instance()->GetComms();
      OpenFile(gmyFileName);
      Geometry gmy = ReadPreamble();
      voxelSize = gmy.GetVoxelSize();
      origin = gmy.GetOrigin();

      HEMELB_MPI_CALL(MPI_File_close, (&file));
      HEMELB_MPI_CALL(MPI_Info_free, (&fileInfo));
    }
  }
}
