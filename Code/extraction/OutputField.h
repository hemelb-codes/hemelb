
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_EXTRACTION_OUTPUTFIELD_H
#define HEMELB_EXTRACTION_OUTPUTFIELD_H

namespace hemelb
{
  namespace extraction
  {
    struct OutputField
    {
        // #658 Refactor out enum
        enum FieldType
        {
          Pressure,
          Velocity,
          ShearStress,
          VonMisesStress,
          ShearRate,
          StressTensor,
          Traction,
          TangentialProjectionTraction,
          MpiRank
        };

        std::string name;
        FieldType type;
    };
  }
}

#endif /* HEMELB_EXTRACTION_OUTPUTFIELD_H */
