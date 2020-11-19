// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_TESTS_HELPERS_LADDFAIL_H
#define HEMELB_TESTS_HELPERS_LADDFAIL_H

#define HEMELB_TESTS_HELPERS_LADDFAIL_Q(str) #str
#define HEMELB_TESTS_HELPERS_LADDFAIL_QQ(str) HEMELB_TESTS_HELPERS_LADDFAIL_Q(str)

#define LADD_FAIL() if (std::strcmp(HEMELB_TESTS_HELPERS_LADDFAIL_QQ(HEMELB_INLET_BOUNDARY), "LADDIOLET") == 0 || std::strcmp(HEMELB_TESTS_HELPERS_LADDFAIL_QQ(HEMELB_OUTLET_BOUNDARY), "LADDIOLET") == 0) FAIL("This test has been marked as fatal when HemeLB is built with LADDIOLET. See ticket #605.")


#endif // ONCE
