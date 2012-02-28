#ifndef HEMELB_EXTRACTION_PROPERTYOUTPUTFILE_H
#define HEMELB_EXTRACTION_PROPERTYOUTPUTFILE_H

#include <string>
#include "extraction/PropertyOutput.h"

namespace hemelb
{
  namespace extraction
  {
    struct PropertyOutputFile
    {
        std::string filename;
        std::vector<PropertyOutput> propertyOutput;
    };
  }
}

#endif /* HEMELB_EXTRACTION_PROPERTYOUTPUTFILE_H */
