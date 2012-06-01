#ifndef HEMELB_GEOMETRY_NEEDS_NEEDS_HPP
#define HEMELB_GEOMETRY_NEEDS_NEEDS_HPP
#include "geometry/needs/Needs.h"
#include "log/Logger.h"
#include "util/utilityFunctions.h"

namespace hemelb
{
  namespace geometry
  {
    template<class Net> NeedsBase<Net>::NeedsBase(const site_t blockCount,
                                                  const std::vector<bool>& readBlock,
                                                  const proc_t areadingGroupSize,
                                                  Net & anet,
                                                  bool shouldValidate) :
        procsWantingBlocksBuffer(blockCount), net(anet), communicator(anet.GetCommunicator()), readingGroupSize(areadingGroupSize), shouldValidate(shouldValidate)
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
      std::vector<int> blocksNeededSizes(communicator.GetSize());

      for (proc_t readingCore = 0; readingCore < readingGroupSize; readingCore++)
      {
        blocksNeededSize[readingCore] = blocksNeededHere[readingCore].size();
        anet.RequestGatherSend(blocksNeededSize[readingCore], readingCore);

      }
      if (communicator.GetRank() < readingGroupSize)
      {
        anet.RequestGatherReceive(blocksNeededSizes);
      }
      anet.Dispatch();
      // Communicate the arrays of needed blocks
      std::vector<std::vector<site_t> > blocksNeededOn(communicator.GetSize());

      for (proc_t readingCore = 0; readingCore < readingGroupSize; readingCore++)
      {

        if (readingCore == communicator.GetRank())
        {
          for (proc_t sendingCore = 0; sendingCore < communicator.GetSize(); sendingCore++)
          {
            blocksNeededOn[sendingCore].resize(blocksNeededSizes[sendingCore]);
          }

        }
        anet.RequestGatherVSend(blocksNeededHere[readingCore], readingCore);
      }
      if (communicator.GetRank() < readingGroupSize)
      {
        anet.RequestGatherVReceive(blocksNeededOn);
      }
      anet.Dispatch();
      if (communicator.GetRank() < readingGroupSize)
      {
        // Transpose the blocks needed on cores matrix
        for (proc_t sendingCore = 0; sendingCore < communicator.GetSize(); sendingCore++)
        {
          for (std::vector<site_t>::iterator need = blocksNeededOn[sendingCore].begin();
              need != blocksNeededOn[sendingCore].end(); need++)
          {
            procsWantingBlocksBuffer[*need].push_back(sendingCore);
          } //for need
        } //for sendingCore
      } //if a reading core

      // If in debug mode, also reduce the old-fashioned way, and check these are the same.
      if (shouldValidate)
      {
        Validate(blockCount, readBlock);
      }
    }

    template<class Net> void NeedsBase<Net>::Validate(const site_t blockCount, const std::vector<bool>& readBlock)
    {
      std::vector<int> procsWantingThisBlockBuffer(communicator.GetSize());
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
                   communicator.GetCommunicator());

        if (communicator.GetRank() == readingCore)
        {

          for (proc_t needingProcOld = 0; needingProcOld < communicator.GetSize(); needingProcOld++)
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
      } // for blocks
    } //constructor

    template<class Net> proc_t NeedsBase<Net>::GetReadingCoreForBlock(const site_t blockNumber) const
    {
      return proc_t(blockNumber % readingGroupSize);
    }
  } //namespace
} //namespace
#endif //HEMELB_GEOMETRY_NEEDS_NEEDS_HPP
