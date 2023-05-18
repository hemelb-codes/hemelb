// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "reporting/BuildInfo.h"
#include "build_info.h"

namespace hemelb::reporting {

    void BuildInfo::Report(Dict &dictionary) {
        Dict build = dictionary.AddSectionDictionary("BUILD");
        build.SetValue("REVISION", build_info::REVISION_HASH);
        build.SetValue("TYPE", build_info::BUILD_TYPE);
        build.SetValue("OPTIMISATION", build_info::OPTIMISATION);
        build.SetBoolValue("USE_SSE3", build_info::USE_SSE3);
        build.SetValue("TIME", build_info::BUILD_TIME);
        build.SetValue("READING_GROUP_SIZE", build_info::READING_GROUP_SIZE);
        build.SetValue("LATTICE_TYPE", build_info::LATTICE);
        build.SetValue("KERNEL_TYPE", build_info::KERNEL);
        build.SetValue("WALL_BOUNDARY_CONDITION", build_info::WALL_BOUNDARY);
        build.SetValue("INLET_BOUNDARY_CONDITION", build_info::INLET_BOUNDARY);
        build.SetValue("OUTLET_BOUNDARY_CONDITION", build_info::OUTLET_BOUNDARY);
        build.SetBoolValue("SEPARATE_CONCERNS", build_info::SEPARATE_CONCERNS);
        build.SetValue("ALLTOALL_IMPLEMENTATION", build_info::ALLTOALL_IMPLEMENTATION);
        build.SetValue("GATHERS_IMPLEMENTATION", build_info::GATHERS_IMPLEMENTATION);
        build.SetValue("POINTPOINT_IMPLEMENTATION", build_info::POINTPOINT_IMPLEMENTATION);
        build.SetValue("STENCIL", build_info::STENCIL);
    }
}