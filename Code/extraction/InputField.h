
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_EXTRACTION_INPUTFIELD_H
#define HEMELB_EXTRACTION_INPUTFIELD_H

namespace hemelb
{
  namespace extraction
  {
    struct InputField
    {
      std::string name;
      uint32_t numberOfFloats;
      double offset;
    };
  }
}

#endif /* HEMELB_EXTRACTION_INPUTFIELD_H */
