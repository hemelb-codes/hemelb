# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

add_library(hemelb_unittests main.cc)
hemelb_add_target_dependency_cppunit(hemelb_unittests)
# ReporterTests directly use ctemplate
hemelb_add_target_dependency_ctemplate(hemelb_unittests)
