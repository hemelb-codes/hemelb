// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <vector>
#include "redblood/WallCellPairIterator.h"
#include "log/Logger.h"

namespace hemelb
{
  namespace redblood
  {
    bool WallCellPairIterator::operator++()
    {
      do
      {
        if (firstWallNode == wallNodes.end())
        {
          return false;
        }
        if (firstCellNode != lastCellNode)
        {
          ++firstCellNode;
        }
        else if (++box_iterator)
        {
          std::tie(firstCellNode, lastCellNode) = cellNodes(firstWallNode->second.node
              + *box_iterator);
        }
        else if (++firstWallNode != wallNodes.end())
        {
          box_iterator = BorderBoxIterator(firstWallNode->second.nearBorder);
          // should always at least include cell of firstWallNode
          assert(static_cast<bool>(box_iterator));
          std::tie(firstCellNode, lastCellNode) = cellNodes(firstWallNode->second.node
              + *box_iterator);
        }
        else
        {
          // reached end of loop
          return false;
        }
      }
      while (not isValid());
      return true;
    }

    bool WallCellPairIterator::isValid() const
    {
      return firstCellNode != lastCellNode and box_iterator
          and (*firstCellNode - firstWallNode->second.node).GetMagnitude() < cutoff;
    }

    std::tuple<DivideConquerCells const&, DivideConquer<hemelb::redblood::WallNode> const&,
        hemelb::LatticeDistance> iterate(
        hemelb::redblood::DivideConquerCells const& cellDnC,
        hemelb::redblood::DivideConquer<hemelb::redblood::WallNode> const& wallDnC,
        hemelb::LatticeDistance const & cutoff)
    {
      using namespace hemelb::redblood;
      return std::make_tuple(std::cref(cellDnC), std::cref(wallDnC), cutoff);
    }

    WallCellPairIterator::WallCellPairIterator(
        DivideConquerCells const& cellNodes, DivideConquer<WallNode> const &wallNodes,
        LatticeDistance cutoff, DivideConquerCells::const_iterator const &firstCellNode,
        DivideConquerCells::const_iterator const &lastCellNode,
        DivideConquer<WallNode>::const_iterator const &firstWallNode) :
        cellNodes(cellNodes), wallNodes(wallNodes), firstCellNode(firstCellNode),
            lastCellNode(lastCellNode), firstWallNode(firstWallNode), box_iterator(0),
            cutoff(cutoff)
    {
      if (firstWallNode == wallNodes.end())
      {
        return;
      }
      try
      {
        box_iterator = BorderBoxIterator(firstWallNode->second.nearBorder);
        assert(static_cast<bool>(box_iterator));
        if (firstCellNode == lastCellNode)
        {
          std::tie(this->firstCellNode, this->lastCellNode) = cellNodes(firstWallNode->second.node
              + *box_iterator);
        }
        if (not isValid())
        {
          operator++();
        }
      }
      catch (...)
      {
        this->firstWallNode = wallNodes.end();
        throw;
      }
    }
  }
} // namespace hemelb::redblood

namespace std
{
  //! Overload so we can work with for-range loop
  hemelb::redblood::WallCellPairIterator begin(
      tuple<hemelb::redblood::DivideConquerCells const&,
          hemelb::redblood::DivideConquer<hemelb::redblood::WallNode> const&,
          hemelb::LatticeDistance> args)
  {
    using namespace hemelb::redblood;
    return
    { get<0>(args), get<1>(args), get<2>(args), WallCellPairIterator::Begin()};
  }

  hemelb::redblood::WallCellPairIterator end(
      tuple<hemelb::redblood::DivideConquerCells const&,
          hemelb::redblood::DivideConquer<hemelb::redblood::WallNode> const&,
          hemelb::LatticeDistance> args)
  {
    using namespace hemelb::redblood;
    return
    { get<0>(args), get<1>(args), get<2>(args), WallCellPairIterator::End()};
  }
} // std
