//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_REDBLOOD_WALL_CELL_PAIR_ITERATOR_H
#define HEMELB_REDBLOOD_WALL_CELL_PAIR_ITERATOR_H

#include <assert.h>
#include "units.h"
#include "redblood/DivideConquer.h"
#include "redblood/CellCell.h"
#include "geometry/LatticeData.h"

namespace hemelb
{
  namespace redblood
  {
    //! References a node of a mesh in the divide-and-conquer box
    class WallNode
    {
      public:
        //! Index of node in mesh
        LatticePosition node;
        //! Whether the node is near the border of the cube;
        size_t nearBorder;
    };

    //! \brief Creates divide and conquer box of wall nodes
    //! \brief The items of the divide and conquer box contain the location of the wall-node (not
    //! the location of the associated fluid site) as well as how close the node is to neighboring
    //! boxes.
    template<class LATTICE>
      DivideConquer<WallNode> createWallNodeDnC(
        geometry::LatticeData const&latticeData, LatticeDistance boxSize,
        LatticeDistance interactionDistance)
      {
        DivideConquer<WallNode> result(boxSize);
        for(site_t i(0); i < latticeData.GetLocalFluidSiteCount(); ++i)
        {
          auto const site = latticeData.GetSite(i);
          if(not site.IsWall())
          {
            continue;
          }
          for (Direction i(1); i < LATTICE::NUMVECTORS; ++i)
          {
            if (not site.HasWall(i))
            {
              continue;
            }
            LatticeDistance const distance = site.GetWallDistance < LATTICE > (i);

            // Direction of streaming from wall to this site
            LatticePosition const direction(LATTICE::CX[i], LATTICE::CY[i], LATTICE::CZ[i]);
            auto const wallnode
              = LatticePosition(site.GetGlobalSiteCoords()) + direction.GetNormalised() * distance;
            auto const nearness = figureNearness(result, wallnode, interactionDistance);
            result.insert(wallnode, {wallnode, nearness});
          }
        }
        return result;
      }
  }
}
#endif
