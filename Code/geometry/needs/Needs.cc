// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#include "geometry/needs/Needs.h"
#include "log/Logger.h"
#include "util/utilityFunctions.h"

namespace hemelb
{
  namespace geometry
  {
    Needs::Needs(const site_t blockCount,
                 const std::vector<bool>& readBlock,
                 const proc_t readingGroupSize,
                 net::InterfaceDelegationNet & net,
                 bool shouldValidate) :
        procsWantingBlocksBuffer(blockCount), net(net), communicator(net.GetCommunicator()), readingGroupSize(readingGroupSize), shouldValidate(shouldValidate)
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
      std::vector<int> blocksNeededSizes(communicator.Size());

      for (proc_t readingCore = 0; readingCore < readingGroupSize; readingCore++)
      {
        blocksNeededSize[readingCore] = blocksNeededHere[readingCore].size();
        net.RequestGatherSend(blocksNeededSize[readingCore], readingCore);

      }
      if (communicator.Rank() < readingGroupSize)
      {
        net.RequestGatherReceive(blocksNeededSizes);
      }
      net.Dispatch();
      // Communicate the arrays of needed blocks

      for (proc_t readingCore = 0; readingCore < readingGroupSize; readingCore++)
      {
        net.RequestGatherVSend(blocksNeededHere[readingCore], readingCore);
      }

      std::vector<site_t> blocksNeededOn;

      if (communicator.Rank() < readingGroupSize)
      {
        net.RequestGatherVReceive(blocksNeededOn, blocksNeededSizes);
      }
      net.Dispatch();
      if (communicator.Rank() < readingGroupSize)
      {
        int needsPassed = 0;
        // Transpose the blocks needed on cores matrix
        for (proc_t sendingCore = 0; sendingCore < communicator.Size(); sendingCore++)
        {
          for (int needForThisSendingCore = 0; needForThisSendingCore < blocksNeededSizes[sendingCore];
              ++needForThisSendingCore)
          {
            procsWantingBlocksBuffer[blocksNeededOn[needsPassed]].push_back(sendingCore);

            ++needsPassed;
          }
        } //for sendingCore
      } //if a reading core

      // If in debug mode, also reduce the old-fashioned way, and check these are the same.
      if (shouldValidate)
      {
        Validate(blockCount, readBlock);
      }
    }

    void Needs::Validate(const site_t blockCount, const std::vector<bool>& readBlock)
    {
      std::vector<int> procsWantingThisBlockBuffer(communicator.Size());
      for (site_t block = 0; block < blockCount; ++block)
      {
        int neededHere = readBlock[block];
        proc_t readingCore = GetReadingCoreForBlock(block);
        MPI_Gather(&neededHere,
                   1,
                   net::MpiDataType<int>(),
                   &procsWantingThisBlockBuffer[0],
                   1,
                   net::MpiDataType<int>(),
                   readingCore,
                   communicator);

        if (communicator.Rank() == readingCore)
        {

          for (proc_t needingProcOld = 0; needingProcOld < communicator.Size(); needingProcOld++)
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
              log::Logger::Log<log::Critical, log::OnePerCore>("Was not expecting block %i to be needed on proc %i, but was",
                                                            block,
                                                            needingProcOld);
            }
            if ( (!found) && procsWantingThisBlockBuffer[needingProcOld] && (needingProcOld != readingCore))
            {
              log::Logger::Log<log::Critical, log::OnePerCore>("Was expecting block %i to be needed on proc %i, but was not",
                                                            block,
                                                            needingProcOld);
            } //if problem
          } // for old proc
        } // if reading group
      } // for blocks
    } //constructor

    proc_t Needs::GetReadingCoreForBlock(const site_t blockNumber) const
    {
      return proc_t(blockNumber % readingGroupSize);
    }
  } //namespace
} //namespace
