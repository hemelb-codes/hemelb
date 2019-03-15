
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "geometry/LatticeData.h"



using namespace hemelb;
using namespace hemelb::geometry;



__global__
void CopyReceivedKernel(
  const site_t* streamingIndicesForReceivedDistributions,
  const distribn_t* fOldShared,
  distribn_t* fNew,
  site_t totalSharedFs
)
{
  int i = blockIdx.x * blockDim.x + threadIdx.x;

  if ( i >= totalSharedFs )
  {
    return;
  }

  fNew[streamingIndicesForReceivedDistributions[i]] = fOldShared[i];
}



void LatticeData::CopyReceivedGPU(int blockSize)
{
  if ( totalSharedFs == 0 )
  {
    return;
  }

  const int GRID_SIZE = (totalSharedFs + blockSize - 1) / blockSize;

  CopyReceivedKernel<<<GRID_SIZE, blockSize>>>(
    streamingIndicesForReceivedDistributions_dev,
    GetFOldGPU(neighbouringProcs[0].FirstSharedDistribution),
    GetFNewGPU(0),
    totalSharedFs
  );
  CUDA_SAFE_CALL(cudaGetLastError());
}
