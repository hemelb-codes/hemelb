
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "units.h"



namespace hemelb {
namespace geometry {



__global__ void LatticeData_CopyReceived(
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



__host__ void LatticeData_CopyReceivedGPU(
  const site_t* streamingIndicesForReceivedDistributions,
  const distribn_t* fOldShared,
  distribn_t* fNew,
  site_t totalSharedFs
)
{
  const int BLOCK_SIZE = 256;
  const int GRID_SIZE = (totalSharedFs + BLOCK_SIZE - 1) / BLOCK_SIZE;

  LatticeData_CopyReceived<<<GRID_SIZE, BLOCK_SIZE>>>(
    streamingIndicesForReceivedDistributions,
    fOldShared,
    fNew,
    totalSharedFs
  );
}



}
}
