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
#include <tuple>
#include <iterator>
#include "units.h"
#include "redblood/DivideConquer.h"
#include "redblood/CellCell.h"
#include "redblood/Borders.h"
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
        LatticeDistance interactionDistance);

    //! Iterates over pairs of wall-node, cell-node which are within a given distance
    class WallCellPairIterator
    {
      public:
        //! Underlying type of the dereferenced object
        struct value_type
        {
          LatticePosition const &cellNode;
          LatticePosition const &wallNode;
        };
        typedef std::unique_ptr<WallCellPairIterator::value_type> pointer;
        typedef WallCellPairIterator::value_type const& reference;
        typedef std::ptrdiff_t difference_type;
        typedef std::input_iterator_tag iterator_category;
        //! Tag to initialize a begining-of-range iterator
        struct Begin
        {
        };
        //! Tag to initialize an end-of-range iterator
        struct End
        {
        };
        //! Creates an iterator of Wall and Cell nodes
        WallCellPairIterator(
            DivideConquerCells const& cellNodes, DivideConquer<WallNode> const &wallNodes,
            LatticeDistance cutoff, WallCellPairIterator::Begin const &)
          : WallCellPairIterator(
              cellNodes, wallNodes, cutoff, cellNodes.end(), cellNodes.end(), wallNodes.begin())
        {
        }
        WallCellPairIterator(
            DivideConquerCells const& cellNodes, DivideConquer<WallNode> const &wallNodes,
            LatticeDistance cutoff)
          : WallCellPairIterator(cellNodes, wallNodes, cutoff, Begin())
        {
        }
        //! Creates an iterator of Wall and Cell nodes
        WallCellPairIterator(
            DivideConquerCells const& cellNodes, DivideConquer<WallNode> const &wallNodes,
            LatticeDistance cutoff, WallCellPairIterator::End const &)
          : WallCellPairIterator(
              cellNodes, wallNodes, cutoff, cellNodes.end(), cellNodes.end(), wallNodes.end())
        {
        }
        //! Positions of the wall and cell nodes
        value_type operator*() const
        {
          assert(static_cast<bool>(*this));
          return {*firstCellNode, firstWallNode->second.node};
        }
        //! Do not use. Only exists to model an InputIterator
        pointer operator->() const
        {
          assert(static_cast<bool>(*this));
          return pointer{new value_type{*firstCellNode, firstWallNode->second.node}};
        }

        bool operator++();
        WallCellPairIterator operator++(int)
        {
          assert(static_cast<bool>(*this));
          WallCellPairIterator result(
              cellNodes, wallNodes, cutoff, firstCellNode, lastCellNode, firstWallNode);
          operator++();
          return result;
        }

        operator bool() const
        {
          return firstWallNode != wallNodes.end() and isValid();
        }

      protected:
        //! Reference to node vertices via divide and conquer boxes
        DivideConquerCells const& cellNodes;
        //! wall node iterator
        DivideConquer<WallNode> const& wallNodes;
        //! cell node iterator
        DivideConquerCells::const_iterator firstCellNode;
        //! cell node iterator
        DivideConquerCells::const_iterator lastCellNode;
        //! wall node iterator
        DivideConquer<WallNode>::const_iterator firstWallNode;
        //! iterates over difference boxes around current wall node
        BorderBoxIterator box_iterator;
        //! Maximum interaction distance
        LatticeDistance cutoff;

        //! Whether current pair is within range
        bool isValid() const;

        //! Creates an iterator of Wall and Cell nodes
        WallCellPairIterator(
            DivideConquerCells const& cellNodes, DivideConquer<WallNode> const &wallNodes,
            LatticeDistance cutoff,
            DivideConquerCells::const_iterator const &firstCellNode,
            DivideConquerCells::const_iterator const &lastCellNode,
            DivideConquer<WallNode>::const_iterator const &firstWallNode);
    };

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

    //! \brief Creates an object that std::begin and std::end can act on
    //! \details This makes it possible to use ranged for loops as well as iterate in more
    //! tradditional manner:
    //! \code{.cpp}
    //!   auto const last = std::end(iterate(cellDnC, wallDnC, interactionDistance));
    //!   auto first = std::begin(iterate(cellDnC, wallDnC, interactionDistance));
    //!   for(; first != last; ++first)
    //!     std::cout << first->wallNode << " " << first->cellNode << "\n";
    //| \endcode
    //! or
    //! \code{.cpp}
    //!   for(auto const nodes: iterate(cellDnC, wallDnC, interactionDistance))
    //!     std::cout << nodes.wallNode << " " << nodes.cellNode << "\n";
    //| \endcode
    std::tuple<
       DivideConquerCells const&,
       DivideConquer<hemelb::redblood::WallNode> const&,
       hemelb::LatticeDistance
     > iterate(
       hemelb::redblood::DivideConquerCells const& cellDnC,
       hemelb::redblood::DivideConquer<hemelb::redblood::WallNode> const& wallDnC,
       hemelb::LatticeDistance const & cutoff
     );

    //! \brief Computes cell <-> cell interactions and spread to grid
    //! \details Given a partition of the cells' nodes and node <-> node interaction
    //! functional, computes the short-range that can occur between cells that are
    //! too close to one another. The interaction forces are computed and spread to
    //! the lattice.
    template <class STENCIL>
    void addCell2WallInteractions(
        DivideConquerCells const &cellDnC,
        DivideConquer<WallNode> const &wallDnC,
        Node2NodeForce const &functional,
        geometry::LatticeData &latticeData
    );
  }
}

namespace std
{
  //! Overload so we can work with for-range loop
  hemelb::redblood::WallCellPairIterator begin(
    tuple<
      hemelb::redblood::DivideConquerCells const&,
      hemelb::redblood::DivideConquer<hemelb::redblood::WallNode> const&,
      hemelb::LatticeDistance
    > args
  );

  //! Overload so we can work with for-range loop
  hemelb::redblood::WallCellPairIterator end(
    tuple<
      hemelb::redblood::DivideConquerCells const&,
      hemelb::redblood::DivideConquer<hemelb::redblood::WallNode> const&,
      hemelb::LatticeDistance
    > args
  );
} // std


namespace hemelb
{
  namespace redblood
  {
    // implementation must happen after begin and end have been declared
    template <class STENCIL>
    void addCell2WallInteractions(
        DivideConquerCells const &cellDnC,
        DivideConquer<WallNode> const &wallDnC,
        Node2NodeForce const &functional,
        geometry::LatticeData &latticeData
    )
    {
      for(auto const nodes: iterate(cellDnC, wallDnC, functional.cutoff))
      {
        auto const force = functional(nodes.cellNode, nodes.wallNode);
        spreadForce<STENCIL>(nodes.cellNode, latticeData, force);
      }
    }
  }
}

#endif
