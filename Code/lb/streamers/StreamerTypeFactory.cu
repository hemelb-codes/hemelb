#include "lb/lattices/Lattices.h"
#include "lb/streamers/StreamerTypeFactory.cuh"
#include "iostream"

#define CUDA_SAFE_CALL(x)                           \
{                                                   \
    cudaError_t error = x;                          \
    if ( error != cudaSuccess ) {                   \
      const char *name = cudaGetErrorName(error);   \
      const char *str = cudaGetErrorString(error);  \
      std::cerr << "\n"                             \
                << "CUDA Error at " #x "\n"         \
                << name << ": " << str << "\n";     \
      exit(1);                                      \
    }                                               \
}


namespace hemelb {
namespace lb {
namespace streamers {

const int MAX_LATTICE_NUMVECTORS=27;
__constant__ Direction d_NUMVECTORS = MAX_LATTICE_NUMVECTORS;
__constant__ int d_CX[MAX_LATTICE_NUMVECTORS];
__constant__ int d_CY[MAX_LATTICE_NUMVECTORS];
__constant__ int d_CZ[MAX_LATTICE_NUMVECTORS];
__constant__ int* d_discreteVelocityVectors[3] = {d_CX, d_CY, d_CZ};
__constant__ distribn_t d_CXD[MAX_LATTICE_NUMVECTORS];
__constant__ distribn_t d_CYD[MAX_LATTICE_NUMVECTORS];
__constant__ distribn_t d_CZD[MAX_LATTICE_NUMVECTORS];
__constant__ distribn_t d_EQMWEIGHTS[MAX_LATTICE_NUMVECTORS];
__constant__ Direction d_INVERSEDIRECTIONS[MAX_LATTICE_NUMVECTORS];

template <typename LatticeType>
void FillGPUConstantMemory() {
Direction num_vectors = LatticeType::NUMVECTORS;
CUDA_SAFE_CALL(cudaMemcpyToSymbol(d_NUMVECTORS, &num_vectors, sizeof(Direction)));
CUDA_SAFE_CALL(cudaMemcpyToSymbol(d_CX, LatticeType::CX, LatticeType::NUMVECTORS * sizeof(int)));
CUDA_SAFE_CALL(cudaMemcpyToSymbol(d_CY, LatticeType::CY, LatticeType::NUMVECTORS * sizeof(int)));
CUDA_SAFE_CALL(cudaMemcpyToSymbol(d_CZ, LatticeType::CZ, LatticeType::NUMVECTORS * sizeof(int)));
CUDA_SAFE_CALL(cudaMemcpyToSymbol(d_discreteVelocityVectors, LatticeType::discreteVelocityVectors, LatticeType::NUMVECTORS * 3 *sizeof(int)));
CUDA_SAFE_CALL(cudaMemcpyToSymbol(d_CXD, LatticeType::CXD, LatticeType::NUMVECTORS * sizeof(distribn_t)));
CUDA_SAFE_CALL(cudaMemcpyToSymbol(d_CYD, LatticeType::CYD, LatticeType::NUMVECTORS * sizeof(distribn_t)));
CUDA_SAFE_CALL(cudaMemcpyToSymbol(d_CZD, LatticeType::CZD, LatticeType::NUMVECTORS * sizeof(distribn_t)));
CUDA_SAFE_CALL(cudaMemcpyToSymbol(d_EQMWEIGHTS, LatticeType::EQMWEIGHTS, LatticeType::NUMVECTORS * sizeof(distribn_t)));
CUDA_SAFE_CALL(cudaMemcpyToSymbol(d_INVERSEDIRECTIONS, LatticeType::INVERSEDIRECTIONS, LatticeType::NUMVECTORS * sizeof(Direction)));
}

template void FillGPUConstantMemory<hemelb::lb::lattices::D3Q15>();
template void FillGPUConstantMemory<hemelb::lb::lattices::D3Q15i>();
template void FillGPUConstantMemory<hemelb::lb::lattices::D3Q19>();
template void FillGPUConstantMemory<hemelb::lb::lattices::D3Q27>();

// type definitions taken from host code
typedef int64_t site_t;
typedef double distribn_t;

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
  distribn_t f[MAX_LATTICE_NUMVECTORS];
  distribn_t density;
  double3 momentum;
  double3 velocity;
  distribn_t f_eq[MAX_LATTICE_NUMVECTORS];
  distribn_t* f_neq = f_eq;
  distribn_t* f_post = f_eq;

  // copy fOld to local memory
  memcpy(&f[0], &fOld[siteIndex * d_NUMVECTORS], d_NUMVECTORS * sizeof(distribn_t));

  // collider.CalculatePreCollision() (collider = Normal, kernel = LBGK)

  // Lattice::CalculateDensityMomentumFEq()
  density = 0.0;
  momentum.x = 0.0;
  momentum.y = 0.0;
  momentum.z = 0.0;

  for ( int j = 0; j < d_NUMVECTORS; ++j )
  {
    density += f[j];
    momentum.x += d_CXD[j] * f[j];
    momentum.y += d_CYD[j] * f[j];
    momentum.z += d_CZD[j] * f[j];
  }

  velocity.x = momentum.x / density;
  velocity.y = momentum.y / density;
  velocity.z = momentum.z / density;

  const distribn_t density_1 = 1. / density;
  const distribn_t momentumMagnitudeSquared =
      momentum.x * momentum.x
      + momentum.y * momentum.y
      + momentum.z * momentum.z;

  for ( int j = 0; j < d_NUMVECTORS; ++j )
  {
    const distribn_t mom_dot_ei =
        d_CX[j] * momentum.x
        + d_CY[j] * momentum.y
        + d_CZ[j] * momentum.z;

    f_eq[j] = d_EQMWEIGHTS[j]
        * (density
            - (3. / 2.) * momentumMagnitudeSquared * density_1
            + (9. / 2.) * density_1 * mom_dot_ei * mom_dot_ei
            + 3. * mom_dot_ei);
  }

  // LBGK::DoCalculateDensityMomentumFeq()
  for ( int j = 0; j < d_NUMVECTORS; ++j )
  {
    f_neq[j] = f[j] - f_eq[j];
  }

  // collider.Collide()

  // LBGK::DoCollide()
  for ( int j = 0; j < d_NUMVECTORS; ++j )
  {
    f_post[j] = f[j] + f_neq[j] * lbmParams_omega;
  }

  // perform streaming
  for ( int j = 0; j < d_NUMVECTORS; ++j )
  {
    if ( Site_HasWall(wallIntersections[siteIndex], j) )
    {
      int outIndex = siteIndex * d_NUMVECTORS + d_INVERSEDIRECTIONS[j];
      fNew[outIndex] = f_post[j];
    }
    else
    {
      int outIndex = neighbourIndices[siteIndex * d_NUMVECTORS + j];
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

  DoStreamAndCollideKernel <<<GRID_SIZE, BLOCK_SIZE>>>(
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
