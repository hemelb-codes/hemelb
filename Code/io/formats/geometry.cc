// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#include "io/formats/geometry.h"

namespace hemelb
{
  namespace io
  {
    namespace formats
    {
      geometry::geometry() :
          displacements(geometry::NumberOfDisplacements)
      {
        displacements[0] = Displacement(-1, -1, -1);
        displacements[1] = Displacement(-1, -1, 0);
        displacements[2] = Displacement(-1, -1, +1);

        displacements[3] = Displacement(-1, 0, -1);
        displacements[4] = Displacement(-1, 0, 0);
        displacements[5] = Displacement(-1, 0, +1);

        displacements[6] = Displacement(-1, +1, -1);
        displacements[7] = Displacement(-1, +1, 0);
        displacements[8] = Displacement(-1, +1, +1);

        displacements[9] = Displacement(0, -1, -1);
        displacements[10] = Displacement(0, -1, 0);
        displacements[11] = Displacement(0, -1, +1);

        displacements[12] = Displacement(0, 0, -1);
        // displacements[] = Displacement( 0, 0, 0);
        displacements[13] = Displacement(0, 0, +1);

        displacements[14] = Displacement(0, +1, -1);
        displacements[15] = Displacement(0, +1, 0);
        displacements[16] = Displacement(0, +1, +1);

        displacements[17] = Displacement(+1, -1, -1);
        displacements[18] = Displacement(+1, -1, 0);
        displacements[19] = Displacement(+1, -1, +1);

        displacements[20] = Displacement(+1, 0, -1);
        displacements[21] = Displacement(+1, 0, 0);
        displacements[22] = Displacement(+1, 0, +1);

        displacements[23] = Displacement(+1, +1, -1);
        displacements[24] = Displacement(+1, +1, 0);
        displacements[25] = Displacement(+1, +1, +1);
      }

      geometry::~geometry()
      {
        // D'tor is trivial here.
      }

      geometry* geometry::singleton = nullptr;
    }
  }
}
