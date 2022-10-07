// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "reporting/BuildInfo.h"
namespace hemelb::reporting {

    void BuildInfo::Report(Dict &dictionary) {
        Dict build = dictionary.AddSectionDictionary("BUILD");
        build.SetValue("REVISION", mercurial_revision_number);
        build.SetValue("TYPE", build_type);
        build.SetValue("OPTIMISATION", optimisation);
        build.SetValue("USE_SSE3", use_sse3);
        build.SetValue("TIME", build_time);
        build.SetValue("READING_GROUP_SIZE", reading_group_size);
        build.SetValue("LATTICE_TYPE", lattice_type);
        build.SetValue("KERNEL_TYPE", kernel_type);
        build.SetValue("WALL_BOUNDARY_CONDITION", wall_boundary_condition);
        build.SetValue("INLET_BOUNDARY_CONDITION", inlet_boundary_condition);
        build.SetValue("OUTLET_BOUNDARY_CONDITION", outlet_boundary_condition);
        build.SetValue("WALL_INLET_BOUNDARY_CONDITION", wall_inlet_boundary_condition);
        build.SetValue("WALL_OUTLET_BOUNDARY_CONDITION", wall_outlet_boundary_condition);
        build.SetValue("SEPARATE_CONCERNS", separate_concerns);
        build.SetValue("ALLTOALL_IMPLEMENTATION", alltoall_impl);
        build.SetValue("GATHERS_IMPLEMENTATION", gathers_impl);
        build.SetValue("POINTPOINT_IMPLEMENTATION", point_point_impl);
        build.SetValue("STENCIL", stencil);
    }
}