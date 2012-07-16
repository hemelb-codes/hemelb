// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_EXTRACTION_PROPERTYOUTPUTFILE_H
#define HEMELB_EXTRACTION_PROPERTYOUTPUTFILE_H

#include <string>
#include <vector>
#include "extraction/GeometrySelector.h"
#include "extraction/OutputField.h"

namespace hemelb
{
  namespace extraction
  {
    struct PropertyOutputFile
    {
        PropertyOutputFile()
        {
          geometry = NULL;
        }

        ~PropertyOutputFile()
        {
          delete geometry;
        }

        std::string filename;
        unsigned long frequency;
        GeometrySelector* geometry;
        std::vector<OutputField> fields;
    };
  }
}

#endif /* HEMELB_EXTRACTION_PROPERTYOUTPUTFILE_H */
