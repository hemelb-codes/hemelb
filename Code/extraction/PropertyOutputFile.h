
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

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
