// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_GEOMETRY_BLOCKTRAVERSERWITHVISITEDBLOCKTRACKER_H
#define HEMELB_GEOMETRY_BLOCKTRAVERSERWITHVISITEDBLOCKTRACKER_H

#include <vector>
#include "geometry/BlockTraverser.h"

namespace hemelb
{
  namespace geometry
  {
    /**
     * On top of BlockTraverser, this class also contains 
     * a record of which blocks have been visited, which
     * is neccessary for the algoritm which uses this. No locations are automatically
     * marked visited, and methods have been created to assist with random access
     * of the lattice data as required by the algorithm
     */
    class BlockTraverserWithVisitedBlockTracker : public geometry::BlockTraverser
    {
      public:
        BlockTraverserWithVisitedBlockTracker(const geometry::LatticeData& iLatDat);

        virtual ~BlockTraverserWithVisitedBlockTracker();

        //Tranverses the block until the next unvisited block is reached.
        //Returns false if the end of the Volume is reached
        bool GoToNextUnvisitedBlock();

        bool IsCurrentBlockVisited() const;

        bool IsBlockVisited(size_t n) const;
        bool IsBlockVisited(util::Vector3D<site_t> n) const;

        void MarkCurrentBlockVisited();

        void MarkBlockVisited(size_t n);
        void MarkBlockVisited(util::Vector3D<site_t> location);

      private:
        std::vector<bool> mBlockVisited;
    };
  }
}

#endif // HEMELB_GEOMETRY_BLOCKTRAVERSERWITHVISITEDBLOCKTRACKER_H
