// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_GEOMETRY_NEEDS_NEEDS_H
#define HEMELB_GEOMETRY_NEEDS_NEEDS_H
#include <vector>
#include "net/net.h"
#include "net/NetworkTopology.h"
namespace hemelb
{
  namespace geometry
  {
    /*
     JH Note: Class created to make the new needs communication strategy in #142 testable.
     */
    /***
     *  Class defining HemeLB needs communication
     Used by geometry reader to know where to send which blocks.
     */
    class Needs
    {
      public:
        /***
         * Constructor for Needs manager.
         * @param BlockCount Count of blocks
         * @param readBlock Which cores need which blocks, as an array of booleans.
         * @param readingGroupSize Number sof cores to use for reading blocks
         * @param net Instance of Net communication class to use.
         */
       Needs(const site_t blockCount,
                          const std::vector<bool>& readBlock,
                          const proc_t readingGroupSize,
                          net::InterfaceDelegationNet &net,
                          bool shouldValidate); // Temporarily during the refactor, constructed just to abstract the block sharing bit

        /***
         * Which processors need a given block?
         * @param block Block number to query
         * @return Vector of ranks in the decomposition topology which need this block
         */
        const std::vector<proc_t> & ProcessorsNeedingBlock(const site_t &block) const
        {
          return procsWantingBlocksBuffer[block];
        }

        /***
         * Which core should be responsible for reading a given block? This core does not necessarily
         * require information about the block
         *
         * @param blockNumber Block number to query
         * @return Rank in the decomposition topology, for core which should read the block.
         */
        proc_t GetReadingCoreForBlock(const site_t blockNumber) const;
      private:
        std::vector<std::vector<proc_t> > procsWantingBlocksBuffer;
        net::InterfaceDelegationNet &net;
        const net::MpiCommunicator & communicator;
        const proc_t readingGroupSize;
        bool shouldValidate;
        void Validate(const site_t blockCount, const std::vector<bool>& readBlock);
    };
  }
}
#endif // HEMELB_GEOMETRY_NEEDS_NEEDS_H
