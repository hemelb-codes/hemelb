// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_REPORTING_BUILDINFO_H
#define HEMELB_REPORTING_BUILDINFO_H
#include "reporting/Reportable.h"

namespace hemelb::reporting
{
    class BuildInfo : public Reportable {
    public:
        void Report(Dict& dictionary) override;
    };
}
#endif
