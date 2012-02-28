#ifndef HEMELB_EXTRACTION_PROPERTYOUTPUT_H
#define HEMELB_EXTRACTION_PROPERTYOUTPUT_H

#include <vector>
#include "extraction/GeometrySelector.h"
#include "extraction/OutputField.h"

namespace hemelb
{
  namespace extraction
  {
    struct PropertyOutput
    {
        ~PropertyOutput()
        {
          delete geometry;
        }

        GeometrySelector* geometry;
        std::vector<OutputField> fields;
        unsigned long frequency;
    };
  }
}

#endif /* HEMELB_EXTRACTION_PROPERTYOUTPUT_H */
