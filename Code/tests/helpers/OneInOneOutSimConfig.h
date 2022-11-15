// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_TESTS_HELPERS_ONEINONEOUTSIMCONFIG_H
#define HEMELB_TESTS_HELPERS_ONEINONEOUTSIMCONFIG_H

#include <string>

#include "io/xml.h"
#include "configuration/SimConfig.h"
#include "configuration/SimBuilder.h"

namespace hemelb::tests::helpers
{
    // TODO: Figure out what this is supposed to be.
    class OneInOneOutSimConfig : public configuration::SimConfig
    {
    public:
        OneInOneOutSimConfig();
    protected:
        void CheckIoletMatchesCMake(const io::xml::Element& ioletEl,
                                    const std::string& requiredBC) const override;
    };

}

#endif
