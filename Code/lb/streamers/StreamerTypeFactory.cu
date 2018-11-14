
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "lb/iolets/InOutLetCosine.cuh"
#include "lb/lattices/D3Q15.cuh"



namespace hemelb {
namespace lb {
namespace streamers {



#define D3Q15 lattices::D3Q15



// geometry/SiteType.h
enum site_type_t
{
  SOLID_TYPE = 0U,
  FLUID_TYPE = 1U,
  INLET_TYPE = 2U,
  OUTLET_TYPE = 3U
};



// geometry/SiteDataBare.h
typedef struct
{
  unsigned wallIntersection;
  unsigned ioletIntersection;
  site_type_t type;
  int ioletId;
} site_data_t;



__device__ bool Site_HasIolet(unsigned ioletIntersection, int direction)
{
  unsigned mask = 1U << (direction - 1);
  return ((ioletIntersection & mask) != 0) && (direction > 0);
}



__device__ bool Site_HasWall(unsigned wallIntersection, int direction)
{
  unsigned mask = 1U << (direction - 1);
  return ((wallIntersection & mask) != 0) && (direction > 0);
}



// lb/lattices/Lattice.h
__device__ void Lattice_CalculateFeq(const distribn_t& density, const double3& momentum, distribn_t* f_eq)
{
  const distribn_t density_1 = 1. / density;
  const distribn_t momentumMagnitudeSquared =
      momentum.x * momentum.x
      + momentum.y * momentum.y
      + momentum.z * momentum.z;

  for ( int j = 0; j < D3Q15::NUMVECTORS; ++j )
  {
    const distribn_t mom_dot_ei =
        D3Q15::CXD[j] * momentum.x
        + D3Q15::CYD[j] * momentum.y
        + D3Q15::CZD[j] * momentum.z;

    f_eq[j] = D3Q15::EQMWEIGHTS[j]
        * (density
            - (3. / 2.) * momentumMagnitudeSquared * density_1
            + (9. / 2.) * density_1 * mom_dot_ei * mom_dot_ei
            + 3. * mom_dot_ei);
  }
}



__global__ void DoStreamAndCollideKernel(
  site_t firstIndex,
  site_t siteCount,
  distribn_t lbmParams_tau,
  distribn_t lbmParams_omega,
  const iolets::InOutLetCosineGPU* inlets,
  const iolets::InOutLetCosineGPU* outlets,
  const site_t* neighbourIndices,
  const site_data_t* siteData,
  const distribn_t* fOld,
  distribn_t* fNew,
  unsigned long timeStep
)
{
  int i = blockIdx.x * blockDim.x + threadIdx.x;

  if ( i >= siteCount )
  {
    return;
  }

  site_t siteIndex = firstIndex + i;

  // initialize hydroVars
  distribn_t f[D3Q15::NUMVECTORS];
  distribn_t density;
  double3 momentum;
  double3 velocity;
  distribn_t f_eq[D3Q15::NUMVECTORS];
  distribn_t* f_neq = f_eq;
  distribn_t* f_post = f_eq;

  // copy fOld to local memory
  memcpy(&f[0], &fOld[siteIndex * D3Q15::NUMVECTORS], D3Q15::NUMVECTORS * sizeof(distribn_t));

  // collider.CalculatePreCollision() (collider = Normal, kernel = LBGK)

  // Lattice::CalculateDensityMomentumFEq()
  density = 0.0;
  momentum.x = 0.0;
  momentum.y = 0.0;
  momentum.z = 0.0;

  for ( int j = 0; j < D3Q15::NUMVECTORS; ++j )
  {
    density += f[j];
    momentum.x += D3Q15::CXD[j] * f[j];
    momentum.y += D3Q15::CYD[j] * f[j];
    momentum.z += D3Q15::CZD[j] * f[j];
  }

  velocity.x = momentum.x / density;
  velocity.y = momentum.y / density;
  velocity.z = momentum.z / density;

  Lattice_CalculateFeq(density, momentum, f_eq);

  // LBGK::DoCalculateDensityMomentumFeq()
  for ( int j = 0; j < D3Q15::NUMVECTORS; ++j )
  {
    f_neq[j] = f[j] - f_eq[j];
  }

  // collider.Collide()

  // LBGK::DoCollide()
  for ( int j = 0; j < D3Q15::NUMVECTORS; ++j )
  {
    f_post[j] = f[j] + f_neq[j] * lbmParams_omega;
  }

  // perform streaming
  site_data_t site = siteData[siteIndex];

  for ( int j = 0; j < D3Q15::NUMVECTORS; ++j )
  {
    if ( Site_HasIolet(site.ioletIntersection, j) )
    {
      // get iolet
      iolets::InOutLetCosineGPU iolet = (site.type == INLET_TYPE)
        ? inlets[site.ioletId]
        : outlets[site.ioletId];

      // get density at the iolet
      distribn_t ghost_density = iolet.GetDensity(timeStep);

      // compute momentum at the iolet
      distribn_t component =
          velocity.x * iolet.normal.x
          + velocity.y * iolet.normal.y
          + velocity.z * iolet.normal.z;

      double3 ghost_momentum;
      ghost_momentum.x = iolet.normal.x * component * ghost_density;
      ghost_momentum.y = iolet.normal.y * component * ghost_density;
      ghost_momentum.z = iolet.normal.z * component * ghost_density;

      // compute f_eq at the iolet
      distribn_t ghost_f_eq[D3Q15::NUMVECTORS];

      Lattice_CalculateFeq(ghost_density, ghost_momentum, ghost_f_eq);

      int outIndex = siteIndex * D3Q15::NUMVECTORS + D3Q15::INVERSEDIRECTIONS[j];
      fNew[outIndex] = ghost_f_eq[D3Q15::INVERSEDIRECTIONS[j]];
    }
    else if ( Site_HasWall(site.wallIntersection, j) )
    {
      int outIndex = siteIndex * D3Q15::NUMVECTORS + D3Q15::INVERSEDIRECTIONS[j];
      fNew[outIndex] = f_post[j];
    }
    else
    {
      int outIndex = neighbourIndices[siteIndex * D3Q15::NUMVECTORS + j];
      fNew[outIndex] = f_post[j];
    }
  }
}



__host__ void DoStreamAndCollideGPU(
  site_t firstIndex,
  site_t siteCount,
  distribn_t lbmParams_tau,
  distribn_t lbmParams_omega,
  const iolets::InOutLetCosineGPU* inlets,
  const iolets::InOutLetCosineGPU* outlets,
  const site_t* neighbourIndices,
  const void* siteData,
  const distribn_t* fOld,
  distribn_t* fNew,
  unsigned long timeStep
)
{
  const int BLOCK_SIZE = 256;
  const int GRID_SIZE = (siteCount + BLOCK_SIZE - 1) / BLOCK_SIZE;

  DoStreamAndCollideKernel<<<GRID_SIZE, BLOCK_SIZE>>>(
    firstIndex,
    siteCount,
    lbmParams_tau,
    lbmParams_omega,
    inlets,
    outlets,
    neighbourIndices,
    (site_data_t*) siteData,
    fOld,
    fNew,
    timeStep
  );
}



}
}
}
