#ifndef HEMELB_GEOMETRY_DECOMPOSITION_HPP
#define HEMELB_GEOMETRY_DECOMPOSITION_HPP
#include "geometry/Decomposition.h"
#include "log/Logger.h"
#include "util/utilityFunctions.h"

namespace hemelb
{
  namespace geometry
  {
    template<class Net> DecompositionBase<Net>::DecompositionBase(const site_t blockCount,
                                                                  const std::vector<bool>& readBlock,
                                                                  const proc_t areadingGroupSize,
                                                                  Net & anet,
                                                                  MPI_Comm comm,
                                                                  const proc_t rank,
                                                                  const proc_t size) :
        procsWantingBlocksBuffer(blockCount), net(anet), decompositionCommunicator(comm), decompositionCommunicatorRank(rank), decompositionCommunicatorSize(size), readingGroupSize(areadingGroupSize)
    {
      // Compile the blocks needed here into an array of indices, instead of an array of bools
      std::vector<std::vector<site_t> > blocksNeededHere(readingGroupSize);
      for (site_t block = 0; block < blockCount; ++block)
      {
        if (readBlock[block])
        {
          blocksNeededHere[GetReadingCoreForBlock(block)].push_back(block);
        }
      }

      // Share the counts of needed blocks
      int blocksNeededSize[readingGroupSize];
      std::vector<int> blocksNeededSizes(decompositionCommunicatorSize);

      for (proc_t readingCore = 0; readingCore < readingGroupSize; readingCore++)
      {
        blocksNeededSize[readingCore] = blocksNeededHere[readingCore].size();

        log::Logger::Log<log::Debug, log::OnePerCore>("Sending count of needed blocks (%i) to core %i from core %i",
                                                      blocksNeededSize[readingCore],
                                                      readingCore,
                                                      decompositionCommunicatorRank);
        anet.RequestGatherSend(blocksNeededSize[readingCore], readingCore);
        if (decompositionCommunicatorRank == readingCore)
        {
          anet.RequestGatherReceive(blocksNeededSizes);
        }

        anet.Dispatch();
      }

      // Communicate the arrays of needed blocks
      std::vector<std::vector<site_t> > blocksNeededOn(decompositionCommunicatorSize);

      for (proc_t readingCore = 0; readingCore < readingGroupSize; readingCore++)
      {

        if (readingCore == decompositionCommunicatorRank)
        {
          for (proc_t sendingCore = 0; sendingCore < decompositionCommunicatorSize; sendingCore++)
          {
            blocksNeededOn[sendingCore].resize(blocksNeededSizes[sendingCore]);
            log::Logger::Log<log::Debug, log::OnePerCore>("Expecting %i needs from core %i",
                                                                    blocksNeededOn[sendingCore].size(),sendingCore);
          }
          anet.RequestGatherVReceive(blocksNeededOn);
          log::Logger::Log<log::Debug, log::OnePerCore>("Receiving lists of blocks needed at core %i",
                                                        decompositionCommunicatorRank);
        }
        anet.RequestGatherVSend(blocksNeededHere[readingCore], readingCore);
        log::Logger::Log<log::Debug, log::OnePerCore>("Sending list of %i needed blocks to core %i from %i",
                                                      blocksNeededHere[readingCore].size(),
                                                      readingCore,
                                                      decompositionCommunicatorRank);

        anet.Dispatch();
      }


      if (decompositionCommunicatorRank < readingGroupSize)
      {
        // Transpose the blocks needed on cores matrix
        for (proc_t sendingCore = 0; sendingCore < decompositionCommunicatorSize; sendingCore++)
        {
          for (std::vector<site_t>::iterator need = blocksNeededOn[sendingCore].begin();
              need != blocksNeededOn[sendingCore].end(); need++)
          {
            procsWantingBlocksBuffer[*need].push_back(sendingCore);
          } //for need
        } //for sendingCore
      } //if a reading core

      // If in debug mode, also reduce the old-fashioned way, and check these are the same.
      if (log::Logger::ShouldDisplay<log::Debug>() && decompositionCommunicator != NULL)
      {

        std::vector<int> procsWantingThisBlockBuffer(decompositionCommunicatorSize);
        for (site_t block = 0; block < blockCount; ++block)
        {
          int neededHere = readBlock[block];
          proc_t readingCore = GetReadingCoreForBlock(block);
          MPI_Gather(&neededHere,
                     1,
                     MpiDataType<int>(),
                     &procsWantingThisBlockBuffer[0],
                     1,
                     MpiDataType<int>(),
                     readingCore,
                     decompositionCommunicator);

          if (decompositionCommunicatorRank == readingCore)
          {

            for (proc_t needingProcOld = 0; needingProcOld < decompositionCommunicatorSize; needingProcOld++)
            {
              bool found = false;
              for (std::vector<proc_t>::iterator needingProc = procsWantingBlocksBuffer[block].begin();
                  needingProc != procsWantingBlocksBuffer[block].end(); needingProc++)
              {
                if (*needingProc == needingProcOld)
                {
                  found = true;
                }
              }
              if (found && (!procsWantingThisBlockBuffer[needingProcOld]))
              {
                log::Logger::Log<log::Debug, log::OnePerCore>("Was not expecting block %i to be needed on proc %i, but was",
                                                              block,
                                                              needingProcOld);
              }
              if ( (!found) && procsWantingThisBlockBuffer[needingProcOld] && (needingProcOld != readingCore))
              {
                log::Logger::Log<log::Debug, log::OnePerCore>("Was expecting block %i to be needed on proc %i, but was not",
                                                              block,
                                                              needingProcOld);
              } //if problem
            } // for old proc
          } // if reading group
        } // for block
      } // if debug mode
    } //constructor

    template<class Net> proc_t DecompositionBase<Net>::GetReadingCoreForBlock(const site_t blockNumber) const
    {
      return proc_t(blockNumber % readingGroupSize);
    }
  } //namespace
} //namespace
#endif //HEMELB_GEOMETRY_DECOMPOSITION_HPP
