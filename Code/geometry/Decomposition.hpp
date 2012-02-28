#ifndef HEMELB_GEOMETRY_DECOMPOSITION_HPP
#define HEMELB_GEOMETRY_DECOMPOSITION_HPP
#include "Decomposition.h"
#include "log/Logger.h"
#include "util/utilityFunctions.h"

namespace hemelb
{
  namespace geometry
  {
    template<class Net> DecompositionBase<Net>::DecompositionBase
    (const site_t blockCount, bool *readBlock, const proc_t readingGroupSize, Net & anet,MPI_Comm comm, const proc_t rank, const proc_t size):
    procsWantingBlocksBuffer(blockCount),net(anet),decompositionCommunicator(comm),decompositionCommunicatorRank(rank),decompositionCommunicatorSize(size)
    {
      // Compile the blocks needed here into an array of indices, instead of an array of bools
      std::vector<site_t> blocks_needed_here;
      for (site_t block = 0; block < blockCount; ++block){
        if (readBlock[block]){
          blocks_needed_here.push_back(block);
        }
      }
      for (proc_t reading_core=0 ; reading_core < readingGroupSize ; reading_core++){
        unsigned int blocks_needed_size=blocks_needed_here.size();
        //if (reading_core==decompositionCommunicatorRank) continue;
        log::Logger::Log<log::Debug, log::OnePerCore>("Sending count of needed blocks (%i) to core %i",blocks_needed_size,reading_core);
        net.RequestSend(&blocks_needed_size, 1, reading_core);
      }

      std::vector<unsigned int>  blocks_needed_sizes(decompositionCommunicatorSize);

      if (decompositionCommunicatorRank < readingGroupSize){
        for (proc_t sending_core=0; sending_core< decompositionCommunicatorSize;sending_core++){
          //if (sending_core==decompositionCommunicatorRank) continue;
          log::Logger::Log<log::Debug, log::OnePerCore>("Receiving count of needed blocks from core %i",sending_core);
          net.RequestReceive(&blocks_needed_sizes[sending_core], 1, sending_core);
        }
      }
      net.Send();
      net.Receive();
      net.Wait();

      // Communicate the needed blocks
      for (proc_t reading_core=0; reading_core<readingGroupSize;reading_core++){
        //if (reading_core==decompositionCommunicatorRank) continue;
        log::Logger::Log<log::Debug, log::OnePerCore>("Sending %i needed blocks to core %i",blocks_needed_here.size(),reading_core);
        net.RequestSend(&blocks_needed_here.front(), blocks_needed_here.size(), reading_core);
      }

      std::vector<std::vector< site_t> >  blocks_needed_on(decompositionCommunicatorSize);

      if (decompositionCommunicatorRank < readingGroupSize){
        for (proc_t sending_core=0; sending_core< decompositionCommunicatorSize;sending_core++){
          blocks_needed_on[sending_core].resize(blocks_needed_sizes[sending_core]);
          //if (sending_core==decompositionCommunicatorRank) continue;
          log::Logger::Log<log::Debug, log::OnePerCore>("Receiving %i needed blocks from core %i",blocks_needed_on[sending_core].size(),sending_core);
          net.RequestReceive(&blocks_needed_on[sending_core].front(), blocks_needed_on[sending_core].size(), sending_core);
        }
      }
      net.Send();
      net.Receive();
      net.Wait();

      if (decompositionCommunicatorRank < readingGroupSize){
        // Transpose the blocks needed on cores matrix
        for (proc_t sending_core=0; sending_core< decompositionCommunicatorSize;sending_core++){
          for (std::vector<site_t>::iterator need=blocks_needed_on[sending_core].begin();need!=blocks_needed_on[sending_core].end();need++){
            for (site_t block = 0; block < blockCount; ++block){
              if (*need==block){
                procsWantingBlocksBuffer[block].push_back(sending_core);
              } //if need the block
            }//for block
          }//for need
        }//for sending_core
      } //if a reading core


      // If in debug mode, also reduce the old-fashioned way, and check these are the same.
      if (log::Logger::ShouldDisplay<log::Debug>())
      {
        
        int* procsWantingThisBlockBuffer = new int[decompositionCommunicatorSize];
        for (site_t block = 0; block < blockCount; ++block)
        {
          int neededHere=readBlock[block];
          proc_t readingCore = proc_t(block % util::NumericalFunctions::min(readingGroupSize, decompositionCommunicatorSize));
          MPI_Gather(&neededHere,1,MpiDataType<int>(),procsWantingThisBlockBuffer,1,MpiDataType<int>(),readingCore,decompositionCommunicator); 
          if (decompositionCommunicatorRank==readingCore)
          {
           
            for (proc_t needing_proc_old=0;needing_proc_old<decompositionCommunicatorSize;needing_proc_old++)
            {
              bool found=false;
              for (std::vector<proc_t>::iterator needing_proc = procsWantingBlocksBuffer[block].begin();
                needing_proc!=procsWantingBlocksBuffer[block].end(); needing_proc++) {
                if (*needing_proc==needing_proc_old) found=true;
              }
              if (found&&(!procsWantingThisBlockBuffer[needing_proc_old])){
                log::Logger::Log<log::Debug, log::OnePerCore>("Was not expecting block %i to be needed on proc %i, but was",
                  block,
                  needing_proc_old
                  );
              }
              if ((!found)&&procsWantingThisBlockBuffer[needing_proc_old])
              {
                log::Logger::Log<log::Debug, log::OnePerCore>("Was expecting block %i to be needed on proc %i, but was not",
                  block,
                  needing_proc_old
                  );
              } //if problem
            } // for old proc
          } // if reading group
        }// for block
        delete[] procsWantingThisBlockBuffer;
      } // if debug mode
    }//constructor
  }//namespace
}//namespace
#endif //HEMELB_GEOMETRY_DECOMPOSITION_HPP
