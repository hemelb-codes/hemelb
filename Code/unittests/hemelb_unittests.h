// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#include "unittests/helpers/helpers.h"
#include "unittests/lbtests/lbtests.h"
#include "unittests/vistests/vistests.h"
#include "unittests/io/io.h"
#include "unittests/reporting/reporting.h"
#include "unittests/configuration/configuration.h"
#include "unittests/geometry/geometry.h"
#include "unittests/SimulationMasterTests.h"
#include "unittests/extraction/extraction.h"
#include "unittests/net/net.h"
#include "unittests/multiscale/multiscale.h"
#ifdef HEMELB_BUILD_MULTISCALE
  #include "unittests/multiscale/mpwide/mpwide.h"
#endif
#include "unittests/redblood.h"
#include "unittests/util/util.h"

#include "unittests/helpers/HasCommsTestFixture.h"
#undef HEMELB_DOING_UNITTESTS
