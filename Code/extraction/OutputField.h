// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_EXTRACTION_OUTPUTFIELD_H
#define HEMELB_EXTRACTION_OUTPUTFIELD_H

namespace hemelb
{
  namespace extraction
  {
    struct OutputField
    {
        enum FieldType
        {
          Pressure,
          Velocity,
          ShearStress,
          VonMisesStress,
          ShearRate
        };

        std::string name;
        FieldType type;
    };
  }
}

#endif /* HEMELB_EXTRACTION_OUTPUTFIELD_H */
