// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
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
#include "unittests/redblood/BendingTests.h"
#include "unittests/redblood/CellTests.h"
#include "unittests/redblood/StencilTests.h"
#include "unittests/redblood/InterpolationTests.h"
#include "unittests/redblood/CellVelocityInterpolTests.h"
#include "unittests/redblood/CellForceSpreadTests.h"
#include "unittests/redblood/Node2NodeTests.h"
#include "unittests/redblood/DivideConquerTests.h"
#include "unittests/redblood/CellArmyTests.h"
#include "unittests/redblood/FlowExtensionTests.h"
#include "unittests/redblood/FaderCellTests.h"
#include "unittests/redblood/VertexBagTests.h"
#include "unittests/redblood/buffer/BufferTests.h"
#include "unittests/redblood/buffer/ColumnsTests.h"
#include "unittests/redblood/CellIOTests.h"
#include "unittests/redblood/CellInserterTests.h"
#include "unittests/redblood/BordersTests.h"
#include "unittests/redblood/WallNodeDnCTests.h"
#include "unittests/redblood/WallCellPairIteratorTests.h"
#include "unittests/redblood/NodeIntegrationTests.h"
#include <unittests/redblood/parallel/NodeParallelizationTests.h>
#include "unittests/redblood/LoadingTimmMeshTests.h"
#include "unittests/redblood/RedBloodMeshVTKDataIOTests.h"
#include "unittests/redblood/LoadDeformedCellTests.h"
#undef HEMELB_DOING_UNITTESTS
#endif
