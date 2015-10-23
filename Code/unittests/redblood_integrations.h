//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//
#ifndef HEMELB_UNITTESTS_REDBLOOD_H
#define HEMELB_UNITTESTS_REDBLOOD_H
namespace hemelb
{
  namespace redblood
  {
    namespace parallel
    {
    }
  }
}
;
namespace hemelb
{
  namespace unittests
  {
    namespace redblood
    {
      using namespace hemelb::redblood;
    }
  }
}
namespace hemelb
{
  namespace unittests
  {
    namespace redblood_parallel
    {
      using namespace hemelb::redblood;
      using namespace hemelb::redblood::parallel;
    }
  }
}
#define HEMELB_DOING_UNITTESTS
#include "unittests/redblood/CellIntegrationTests.h"
#include "unittests/redblood/FadeInOutIntegrationTests.h"
#include "unittests/redblood/ICCBSimulations.h"
#include "unittests/redblood/SadCellTest.h"
#include "unittests/redblood/SadCellIntegrationTests.h"
#include "unittests/redblood/MulticellBenchmarkTests.h"
#include "unittests/redblood/FedosovValidationTests.h"
#undef HEMELB_DOING_UNITTESTS
#endif

