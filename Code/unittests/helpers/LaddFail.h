// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_HELPERS_LADDFAIL_H
#define HEMELB_UNITTESTS_HELPERS_LADDFAIL_H

#include "reporting/BuildInfo.h"
#define HEMELB_UNITTESTS_HELPERS_LADDFAIL_Q(str) #str
#define HEMELB_UNITTESTS_HELPERS_LADDFAIL_QQ(str) HEMELB_UNITTESTS_HELPERS_LADDFAIL_Q(str)

#define LADD_FAIL()                                                     \
  if(::hemelb::reporting::inlet_boundary_condition == "LADDIOLET"       \
      || ::hemelb::reporting::outlet_boundary_condition == "LADDIOLET") \
      CPPUNIT_ASSERT_MESSAGE("This test has been marked as fatal when HemeLB is built with LADDIOLET. See ticket #406.", 1 == 0)

#endif /* HEMELB_UNITTESTS_HELPERS_LADDFAIL_H */
