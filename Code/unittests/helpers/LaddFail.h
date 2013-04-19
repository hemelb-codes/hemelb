//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_HELPERS_LADDFAIL_H
#define HEMELB_UNITTESTS_HELPERS_LADDFAIL_H

#define HEMELB_UNITTESTS_HELPERS_LADDFAIL_Q(str) #str
#define HEMELB_UNITTESTS_HELPERS_LADDFAIL_QQ(str) HEMELB_UNITTESTS_HELPERS_LADDFAIL_Q(str)

#define LADD_FAIL() if (std::strcmp(HEMELB_UNITTESTS_HELPERS_LADDFAIL_QQ(HEMELB_INLET_BOUNDARY), "LADDIOLET") == 0 || std::strcmp(HEMELB_UNITTESTS_HELPERS_LADDFAIL_QQ(HEMELB_OUTLET_BOUNDARY), "LADDIOLET") == 0) CPPUNIT_ASSERT_MESSAGE("This test has been marked as fatal when HemeLB is built with LADDIOLET. See ticket #605.", 1 == 0)


#endif /* HEMELB_UNITTESTS_HELPERS_LADDFAIL_H */
