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
