
// units.h
typedef int64_t site_t;
typedef double distribn_t;



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



// lb/iolets/InOutLetCosine.h
typedef struct
{
  distribn_t minimumSimulationDensity;
  double3 normal;
  double densityMean;
  double densityAmp;
  double phase;
  double period;
  unsigned int warmUpLength;
} iolet_cosine_t;



__device__ distribn_t InOutLetCosine_GetDensity(const iolet_cosine_t& iolet, unsigned long timeStep)
{
  distribn_t w = 2.0 * M_PI / iolet.period;

  distribn_t target = iolet.densityMean + iolet.densityAmp * cos(w * timeStep + iolet.phase);

  if (timeStep >= iolet.warmUpLength)
  {
    return target;
  }

  double interpolationFactor = ((double) timeStep) / ((double) iolet.warmUpLength);

  return interpolationFactor * target + (1. - interpolationFactor) * iolet.minimumSimulationDensity;
}



// lb/lattices/D3Q15.h
namespace D3Q15
{
  __constant__ const int NUMVECTORS = 15;

  __constant__ const distribn_t CXD[] = { 0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1 };
  __constant__ const distribn_t CYD[] = { 0, 0, 0, 1, -1, 0, 0, 1, -1, 1, -1, -1, 1, -1, 1 };
  __constant__ const distribn_t CZD[] = { 0, 0, 0, 0, 0, 1, -1, 1, -1, -1, 1, 1, -1, -1, 1 };

  __constant__ const distribn_t EQMWEIGHTS[] = {
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

  __constant__ const int INVERSEDIRECTIONS[] = { 0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13 };
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
  const iolet_cosine_t* inlets,
  const iolet_cosine_t* outlets,
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
      iolet_cosine_t iolet = (site.type == INLET_TYPE)
        ? inlets[site.ioletId]
        : outlets[site.ioletId];

      // get density at the iolet
      distribn_t ghost_density = InOutLetCosine_GetDensity(iolet, timeStep);

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



namespace hemelb {
namespace lb {
namespace streamers {



__host__ void DoStreamAndCollideGPU(
  site_t firstIndex,
  site_t siteCount,
  distribn_t lbmParams_tau,
  distribn_t lbmParams_omega,
  const iolet_cosine_t* inlets,
  const iolet_cosine_t* outlets,
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
