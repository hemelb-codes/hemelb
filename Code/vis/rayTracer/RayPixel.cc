#include "vis/rayTracer/RayPixel.h"

namespace hemelb
{
  template<>
  MPI_Datatype MpiDataTypeTraits<hemelb::vis::raytracer::RayPixel>::RegisterMpiDataType()
  {
    MPI_Datatype ret = vis::raytracer::RayPixel::GetMPIType();
    MPI_Type_commit(&ret);
    return ret;
  }
}
