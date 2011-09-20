#include "mpiInclude.h"
#include "vis/ColPixel.h"
#include "vis/rayTracer/RayDataNormal.h"
#include "vis/rayTracer/RayDataEnhanced.h"

namespace hemelb
{

  template<>
  MPI_Datatype MpiDataTypeTraits<hemelb::vis::ColPixel<hemelb::vis::raytracer::RayDataNormal> >::RegisterMpiDataType()
  {
    int col_pixel_count = 9;
    int col_pixel_blocklengths[9] = { 1, 1, 1, 1, 1, 1, 1, 1, 1 };

    MPI_Datatype col_pixel_types[9] = { 
      MPI_UNSIGNED,
      MPI_UNSIGNED,
      MPI_UNSIGNED,
      MPI_UNSIGNED,
      MpiDataTypeTraits<hemelb::vis::raytracer::RayDataNormal>
      ::GetMpiDataType(),
      MPI_FLOAT,
      MPI_FLOAT,
      MPI_INT,
      MPI_UB };
    
    MPI_Aint col_pixel_disps[9];

    col_pixel_disps[0] = 0;

    for (int i = 1; i < col_pixel_count; i++)
    {
      if (col_pixel_types[i - 1] == MPI_FLOAT)
      {
        col_pixel_disps[i] = col_pixel_disps[i - 1]
	  + (sizeof(float) * col_pixel_blocklengths[i - 1]);
      }
      else if (col_pixel_types[i - 1] == MPI_INT)
      {
        col_pixel_disps[i] = col_pixel_disps[i - 1] + (sizeof(int) * col_pixel_blocklengths[i - 1]);
      }
      else if (col_pixel_types[i - 1] == MPI_UNSIGNED)
      {
	col_pixel_disps[i] = col_pixel_disps[i - 1] + (sizeof(unsigned) * col_pixel_blocklengths[i - 1]);
      }
      else if (col_pixel_types[i - 1] == MpiDataTypeTraits<hemelb::vis::raytracer::RayDataNormal>
	       ::GetMpiDataType())
      {
	col_pixel_disps[i] = col_pixel_disps[i - 1] + (sizeof(hemelb::vis::raytracer::RayDataNormal)
						       * col_pixel_blocklengths[i - 1]);
      }
      
    }
    MPI_Datatype type;
    MPI_Type_struct(col_pixel_count,
                    col_pixel_blocklengths,
                    col_pixel_disps,
                    col_pixel_types,
                    &type);
    MPI_Type_commit(&type);
    return type;
  }

  template<>
  MPI_Datatype MpiDataTypeTraits<hemelb::vis::ColPixel<hemelb::vis::raytracer::RayDataEnhanced<true> > >::RegisterMpiDataType()
  {
    int col_pixel_count = 9;
    int col_pixel_blocklengths[9] = { 1, 1, 1, 1, 1, 1, 1, 1, 1 };

    MPI_Datatype col_pixel_types[9] = { 
      MPI_UNSIGNED,
      MPI_UNSIGNED,
      MPI_UNSIGNED,
      MPI_UNSIGNED,
      MpiDataTypeTraits<hemelb::vis::raytracer::RayDataEnhanced<true> >
      ::GetMpiDataType(),
      MPI_FLOAT,
      MPI_FLOAT,
      MPI_INT,
      MPI_UB };
    
    MPI_Aint col_pixel_disps[9];

    col_pixel_disps[0] = 0;

    for (int i = 1; i < col_pixel_count; i++)
    {
      if (col_pixel_types[i - 1] == MPI_FLOAT)
      {
        col_pixel_disps[i] = col_pixel_disps[i - 1]
	  + (sizeof(float) * col_pixel_blocklengths[i - 1]);
      }
      else if (col_pixel_types[i - 1] == MPI_INT)
      {
        col_pixel_disps[i] = col_pixel_disps[i - 1] + (sizeof(int) * col_pixel_blocklengths[i - 1]);
      }
      else if (col_pixel_types[i - 1] == MPI_UNSIGNED)
      {
	col_pixel_disps[i] = col_pixel_disps[i - 1] + (sizeof(unsigned) * col_pixel_blocklengths[i - 1]);
      }
      else if (col_pixel_types[i - 1] == MpiDataTypeTraits<hemelb::vis::raytracer::RayDataEnhanced<true> >
	       ::GetMpiDataType())
      {
	col_pixel_disps[i] = col_pixel_disps[i - 1] + (sizeof(hemelb::vis::raytracer::RayDataEnhanced<true>)
						       * col_pixel_blocklengths[i - 1]);
      }
      
    }
    MPI_Datatype type;
    MPI_Type_struct(col_pixel_count,
                    col_pixel_blocklengths,
                    col_pixel_disps,
                    col_pixel_types,
                    &type);
    MPI_Type_commit(&type);
    return type;
  }

template<>
  MPI_Datatype MpiDataTypeTraits<hemelb::vis::ColPixel<hemelb::vis::raytracer::RayDataEnhanced<false> > >::RegisterMpiDataType()
  {
    int col_pixel_count = 9;
    int col_pixel_blocklengths[9] = { 1, 1, 1, 1, 1, 1, 1, 1, 1 };

    MPI_Datatype col_pixel_types[9] = { 
      MPI_UNSIGNED,
      MPI_UNSIGNED,
      MPI_UNSIGNED,
      MPI_UNSIGNED,
      MpiDataTypeTraits<hemelb::vis::raytracer::RayDataEnhanced<false> >
      ::GetMpiDataType(),
      MPI_FLOAT,
      MPI_FLOAT,
      MPI_INT,
      MPI_UB };
    
    MPI_Aint col_pixel_disps[9];

    col_pixel_disps[0] = 0;

    for (int i = 1; i < col_pixel_count; i++)
    {
      if (col_pixel_types[i - 1] == MPI_FLOAT)
      {
        col_pixel_disps[i] = col_pixel_disps[i - 1]
	  + (sizeof(float) * col_pixel_blocklengths[i - 1]);
      }
      else if (col_pixel_types[i - 1] == MPI_INT)
      {
        col_pixel_disps[i] = col_pixel_disps[i - 1] + (sizeof(int) * col_pixel_blocklengths[i - 1]);
      }
      else if (col_pixel_types[i - 1] == MPI_UNSIGNED)
      {
	col_pixel_disps[i] = col_pixel_disps[i - 1] + (sizeof(unsigned) * col_pixel_blocklengths[i - 1]);
      }
      else if (col_pixel_types[i - 1] == MpiDataTypeTraits<hemelb::vis::raytracer::RayDataEnhanced<false> >
	       ::GetMpiDataType())
      {
	col_pixel_disps[i] = col_pixel_disps[i - 1] + (sizeof(hemelb::vis::raytracer::RayDataEnhanced<false>)
						       * col_pixel_blocklengths[i - 1]);
      }
      
    }
    MPI_Datatype type;
    MPI_Type_struct(col_pixel_count,
                    col_pixel_blocklengths,
                    col_pixel_disps,
                    col_pixel_types,
                    &type);
    MPI_Type_commit(&type);
    return type;
  }

}
