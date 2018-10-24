
namespace hemelb {
namespace lb {
namespace streamers {

// type definitions taken from host code
typedef int64_t site_t;
typedef double distribn_t;

// constants (to be implemented as template parameters)
__device__ const int D3Q15_NUMVECTORS = 15;

__device__ const int D3Q15_CX[] = { 0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1 };
__device__ const int D3Q15_CY[] = { 0, 0, 0, 1, -1, 0, 0, 1, -1, 1, -1, -1, 1, -1, 1 };
__device__ const int D3Q15_CZ[] = { 0, 0, 0, 0, 0, 1, -1, 1, -1, -1, 1, 1, -1, -1, 1 };

__device__ const distribn_t D3Q15_CXD[] = { 0.0, 1.0, -1.0, 0.0,  0.0, 0.0,  0.0, 1.0, -1.0,  1.0, -1.0,  1.0, -1.0,  1.0, -1.0};
__device__ const distribn_t D3Q15_CYD[] = { 0.0, 0.0,  0.0, 1.0, -1.0, 0.0,  0.0, 1.0, -1.0,  1.0, -1.0, -1.0,  1.0, -1.0,  1.0};
__device__ const distribn_t D3Q15_CZD[] = { 0.0, 0.0,  0.0, 0.0,  0.0, 1.0, -1.0, 1.0, -1.0, -1.0,  1.0,  1.0, -1.0, -1.0,  1.0};

__device__ const distribn_t D3Q15_EQMWEIGHTS[] = {
  2.0 / 9.0,
  1.0 / 9.0,
  1.0 / 9.0,
  1.0 / 9.0,
  1.0 / 9.0,
  1.0 / 9.0,
  1.0 / 9.0,
  1.0 / 72.0,
  1.0 / 72.0,
  1.0 / 72.0,
  1.0 / 72.0,
  1.0 / 72.0,
  1.0 / 72.0,
  1.0 / 72.0,
  1.0 / 72.0
};

__device__ const int D3Q15_INVERSEDIRECTIONS[] = { 0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13 };



__device__ bool Site_HasWall(unsigned wallIntersection, int direction)
{
  unsigned mask = 1U << max(0, direction - 1);
  return (wallIntersection & mask) != 0;
}



__global__ void DoStreamAndCollideKernel(
  site_t firstIndex,
  site_t siteCount,
  distribn_t lbmParams_tau,
  distribn_t lbmParams_omega,
  const site_t* neighbourIndices,
  const unsigned* wallIntersections,
  const distribn_t* fOld,
  distribn_t* fNew
)
{
  int i = blockIdx.x * blockDim.x + threadIdx.x;

  if ( i >= siteCount )
  {
    return;
  }

  site_t siteIndex = firstIndex + i;

  // initialize hydroVars
  distribn_t f[D3Q15_NUMVECTORS];
  distribn_t density;
  double3 momentum;
  double3 velocity;
  distribn_t f_eq[D3Q15_NUMVECTORS];
  distribn_t* f_neq = f_eq;
  distribn_t* f_post = f_eq;

  // copy fOld to local memory
  memcpy(&f[0], &fOld[siteIndex * D3Q15_NUMVECTORS], D3Q15_NUMVECTORS * sizeof(distribn_t));

  // collider.CalculatePreCollision() (collider = Normal, kernel = LBGK)

  // Lattice::CalculateDensityMomentumFEq()
  density = 0.0;
  momentum.x = 0.0;
  momentum.y = 0.0;
  momentum.z = 0.0;

  for ( int j = 0; j < D3Q15_NUMVECTORS; ++j )
  {
    density += f[j];
    momentum.x += D3Q15_CXD[j] * f[j];
    momentum.y += D3Q15_CYD[j] * f[j];
    momentum.z += D3Q15_CZD[j] * f[j];
  }

  velocity.x = momentum.x / density;
  velocity.y = momentum.y / density;
  velocity.z = momentum.z / density;

  const distribn_t density_1 = 1. / density;
  const distribn_t momentumMagnitudeSquared =
      momentum.x * momentum.x
      + momentum.y * momentum.y
      + momentum.z * momentum.z;

  for ( int j = 0; j < D3Q15_NUMVECTORS; ++j )
  {
    const distribn_t mom_dot_ei =
        D3Q15_CX[j] * momentum.x
        + D3Q15_CY[j] * momentum.y
        + D3Q15_CZ[j] * momentum.z;

    f_eq[j] = D3Q15_EQMWEIGHTS[j]
        * (density
            - (3. / 2.) * momentumMagnitudeSquared * density_1
            + (9. / 2.) * density_1 * mom_dot_ei * mom_dot_ei
            + 3. * mom_dot_ei);
  }

  // LBGK::DoCalculateDensityMomentumFeq()
  for ( int j = 0; j < D3Q15_NUMVECTORS; ++j )
  {
    f_neq[j] = f[j] - f_eq[j];
  }

  // collider.Collide()

  // LBGK::DoCollide()
  for ( int j = 0; j < D3Q15_NUMVECTORS; ++j )
  {
    f_post[j] = f[j] + f_neq[j] * lbmParams_omega;
  }

  // perform streaming
  for ( int j = 0; j < D3Q15_NUMVECTORS; ++j )
  {
    if ( Site_HasWall(wallIntersections[siteIndex], j) )
    {
      int outIndex = siteIndex * D3Q15_NUMVECTORS + D3Q15_INVERSEDIRECTIONS[j];
      fNew[outIndex] = f_post[j];
    }
    else
    {
      int outIndex = neighbourIndices[siteIndex * D3Q15_NUMVECTORS + j];
      fNew[outIndex] = f_post[j];
    }
  }
}



void DoStreamAndCollideGPU(
  site_t firstIndex,
  site_t siteCount,
  distribn_t lbmParams_tau,
  distribn_t lbmParams_omega,
  const site_t* neighbourIndices,
  const unsigned* wallIntersections,
  const distribn_t* fOld,
  distribn_t* fNew
)
{
  const int BLOCK_SIZE = 256;
  const int GRID_SIZE = (siteCount + BLOCK_SIZE - 1) / BLOCK_SIZE;

  DoStreamAndCollideKernel<<<GRID_SIZE, BLOCK_SIZE>>>(
    firstIndex,
    siteCount,
    lbmParams_tau,
    lbmParams_omega,
    neighbourIndices,
    wallIntersections,
    fOld,
    fNew
  );
}



}
}
}
