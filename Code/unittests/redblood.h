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
#include "unittests/redblood/CellCellInteractionTests.h"
#include "unittests/redblood/CellCellInteractionWithGridTests.h"
#include "unittests/redblood/RedBloodMeshTests.h"
#include "unittests/redblood/RedBloodMeshDataIOTests.h"
#include "unittests/redblood/TopologyTests.h"
#include "unittests/redblood/FacetTests.h"
#include "unittests/redblood/EnergyTests.h"
#include "unittests/redblood/GradientTests.h"
#include "unittests/redblood/GradientKernTests.h"
#include "unittests/redblood/CellTests.h"
#include "unittests/redblood/StencilTests.h"
#include "unittests/redblood/InterpolationTests.h"
#include "unittests/redblood/CellVelocityInterpolTests.h"
#include "unittests/redblood/CellForceSpreadTests.h"
#include "unittests/redblood/CellForceSpreadWithWallTests.h"
#include "unittests/redblood/Node2NodeTests.h"
#include "unittests/redblood/DivideConquerTests.h"
#include "unittests/redblood/CellArmyTests.h"
#include "unittests/redblood/CellIntegrationTests.h"
#include "unittests/redblood/FadeInOutIntegrationTests.h"
#include "unittests/redblood/FlowExtensionTests.h"
#include "unittests/redblood/FaderCellTests.h"
#include "unittests/redblood/VertexBagTests.h"
#include "unittests/redblood/buffer/BufferTests.h"
#include "unittests/redblood/buffer/ColumnsTests.h"
#include "unittests/redblood/CellIOTests.h"
#include "unittests/redblood/CellInserterTests.h"
#undef HEMELB_DOING_UNITTESTS
#endif
