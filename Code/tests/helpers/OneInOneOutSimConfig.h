// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_TESTS_HELPERS_ONEINONEOUTSIMCONFIG_H
#define HEMELB_TESTS_HELPERS_ONEINONEOUTSIMCONFIG_H

#include <string>

#include "io/xml/XmlAbstractionLayer.h"
#include "configuration/SimConfig.h"

namespace hemelb
{
  namespace tests
  {
    namespace helpers
    {
      // TODO: Figure out what this is supposed to be.
      class OneInOneOutSimConfig : public configuration::SimConfig
      {
      public:
	OneInOneOutSimConfig(const std::string& path);
      protected:
	virtual void CheckIoletMatchesCMake(const io::xml::Element& ioletEl,
					    const std::string& requiredBC);
      };
    }
  }
}

#endif // ONCE
