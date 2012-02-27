#ifndef HEMELB_GEOMETRY_DECOMPOSITION_HPP
#define HEMELB_GEOMETRY_DECOMPOSITION_HPP
#include "Decomposition.h"
namespace hemelb
{
  namespace geometry
  {
    Decomposition::Decomposition(const site_t blockCount, bool *readBlock, const proc_t readingGroupSize, const topology::NetworkTopology & atopology):
      procsWantingBlocksBuffer(blockCount),net(net::Net(MPI_COMM_WORLD)),topology(atopology)
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
        net.RequestSend(&blocks_needed_size, 1, reading_core);
      }
      std::vector<unsigned int>  blocks_needed_sizes(topology.GetProcessorCount());
      if (topology.GetLocalRank() < readingGroupSize){
        for (proc_t sending_core=0; sending_core< topology.GetProcessorCount();sending_core++){
          net.RequestReceive(&blocks_needed_sizes[sending_core], 1, sending_core);
        }
      }
      net.Send();
      net.Receive();
      net.Wait();

      // Communicate the needed blocks
      for (proc_t reading_core=0; reading_core<readingGroupSize;reading_core++){
        net.RequestSend(&blocks_needed_here.front(), blocks_needed_here.size(), reading_core);
      }
      std::vector<std::vector< site_t> >  blocks_needed_on(topology.GetProcessorCount());
      if (topology.GetLocalRank() < readingGroupSize){
        for (proc_t sending_core=0; sending_core< topology.GetProcessorCount();sending_core++){
          blocks_needed_on[sending_core].resize(blocks_needed_sizes[sending_core]);
          net.RequestReceive(&blocks_needed_on[sending_core].front(), blocks_needed_here.size(), sending_core);
        }
      }
      net.Send();
      net.Receive();
      net.Wait();

      if (topology.GetLocalRank() < readingGroupSize){
        // Transpose the blocks needed on cores matrix
        for (proc_t sending_core=0; sending_core< topology.GetProcessorCount();sending_core++){
          for (std::vector<site_t>::iterator need=blocks_needed_on[sending_core].begin();need!=blocks_needed_on[sending_core].end();need++){
            for (site_t block = 0; block < blockCount; ++block){
              if (*need==block){
                procsWantingBlocksBuffer[block].push_back(sending_core);
              } //if need the block
            }//for block
          }//for need
        }//for sending_core
      } //if a reading core
      
    }//constructor
  }//namespace
}//namespace
#endif //HEMELB_GEOMETRY_DECOMPOSITION_HPP