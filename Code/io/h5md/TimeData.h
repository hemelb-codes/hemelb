#ifndef HEMELB_IO_H5MD_TIMEDATA_H
#define HEMELB_IO_H5MD_TIMEDATA_H

#include <cstdlib>
#include <string>
#include <H5Cpp.h>

namespace hemelb
{
  namespace io
  {
    namespace h5md
    {

      /**
       * H5MD time dependent data.
       *
       * Contains a dataset "value" of type enumeration, string, float or
       * integer (H5::EnumType or subclasses of H5::AtomType).  Rank and
       * dimensionality of the data are fixed at creation time.  The actual rank
       * of the value dataset will be one greater than specified to accommodate
       * a variable dimension allowing the accumulation of samples over time.
       *
       * Contains an integer valued step dataset of rank 1, dimension variable
       * (matches first dimension of value dataset); and a value dataset of type
       * enumeration, string, integer or float of varying rank and dimensions.
       *
       * If time is explicitly store then it will also contain a float-valued
       * dataset of the same rank and dimensionality as step, containing the
       * absolute values of the timesteps sampled.
       */
      class TimeData
      {
        public:

          /**
           * Creates a H5MD Time-Dependent Data group with explicit storage of
           * timesteps.
           *
           * @param location  where to create the group (typically a particles
           *                  or observables group)
           * @param name      the name of the group
           * @param rank      the rank of the data (the values dataset will be
           *                  one greater to accommodate the timestep dimension)
           * @param dims      the dimensions of the data.  May have one or more
           *                  unlimited dimensions (in which case the dataset
           *                  must be chunked).  See HDF5 docs.
           * @param max_dims  maximum dimensions of the data (default the same
           *                  as dims)
           * @param chunk_dims  chunking dimensions (default the same as dims)
           * @param value_type  the HDF5 type of the data (one of enum, integer,
           *                    float or string, default native double)
           * @param step_type   the HDF5 type of the timestep indices (default
           *                    native unsigned long)
           * @param time_type   the HDF5 type of the timestep values (default
           *                    native double)
           * @param timestep_chunk  chunking dimension for timestep values
           *                        (default 1)
           * @param trackTime  whether to create an explicit time dataset to
           *                   store time values or to only store step indices
           *                   to look up values in a shared dataset created
           *                   elsewhere (default omit time dataset)
           * @return a Time-Dependent Data group created at the specified
           *   location.
           */
          static TimeData CreateExplicit(const H5::Group &,
                                         const std::string &,
                                         int, const hsize_t * dims,
                                         const hsize_t * = dims,
                                         const hsize_t * = dims,
                                         const H5::DataType & = H5::PredType::NATIVE_DOUBLE,
                                         const H5::DataType & = H5::PredType::NATIVE_ULONG,
                                         const H5::DataType & = H5::PredType::NATIVE_DOUBLE,
                                         hsize_t = 1,
                                         bool = false);

          /**
           * Creates a H5MD Time-Dependent Data group with implicit storage of
           * timesteps.
           *
           * @param location  where to create the group (typically a particles
           *                  or observables group)
           * @param name      the name of the group
           * @param rank      the rank of the data (the values dataset will be
           *                  one greater to accommodate the timestep dimension)
           * @param dims      the dimensions of the data.  May have one or more
           *                  unlimited dimensions (in which case the dataset
           *                  must be chunked).  See HDF5 docs.
           * @param step_offset  the first timestep (default 0)
           * @param time_offset  the first timestep value (default 0.0)
           * @param max_dims  maximum dimensions of the data (default the same
           *                  as dims)
           * @param chunk_dims  chunking dimensions (default the same as dims)
           * @param value_type  the HDF5 type of the data (one of enum, integer,
           *                    float or string, default native double)
           * @param step_type   the HDF5 type of the timestep indices (default
           *                    native unsigned long)
           * @param time_type   the HDF5 type of the timestep values (default
           *                    native double)
           * @param timestep_chunk  chunking dimension for timestep values
           *                        (default 1)
           * @param trackTime  whether to create an explicit time dataset to
           *                   store time values or to only store step indices
           *                   to look up values in a shared dataset created
           *                   elsewhere (default omit time dataset)
           * @return a Time-Dependent Data group created at the specified
           *   location.
           */
          static TimeData CreateFixed(const H5::Group &,
                                      const std::string &,
                                      int, const hsize_t * dims,
                                      const hsize_t * = dims,
                                      const hsize_t * = dims,
                                      unsigned long = 0, double = 0.0,
                                      const H5::DataType & = H5::PredType::NATIVE_DOUBLE,
                                      const H5::DataType & = H5::PredType::NATIVE_ULONG,
                                      const H5::DataType & = H5::PredType::NATIVE_DOUBLE,
                                      hsize_t = 1,
                                      bool = false);

          static TimeData Open(const H5::Group &, const std::string &);

          bool IsStepExplicit() const;

          bool HasTime() const;

        private:
          TimeData(H5::Group & group, const std::string & name, const H5::DataSet & value,
                   const H5::DataSet & step, const H5::DataSet & time) :
              group(group), name(name), value(value), step(step), time(time)
          {
          }

          /** This group */
          H5::Group group;

          /** The name of this group */
          std::string name;

          /**
           * Contains actual data values from the simulation.
           *
           * Rank is one plus the dimensionality of the data to account for the
           * timestep dimension.  The timestep dimension is unlimited.
           */
          H5::DataSet value;

          /**
           * If timestep storage is explicit then step will have a simple
           * dataspace of rank 1 and dimension variable.  It will have integer
           * type.
           *
           * If timestep storage is fixed then step will have a scalar dataspace
           * and will contain the integer increment between timesteps for values
           * stored.  It may have an attribute, offset, containing the first
           * integer timestep value sampled.
           */
          H5::DataSet step;

          /**
           * If timestep storage is explicit then time will have a simple
           * dataspace of rank 1 and dimension variable.  It will have floating
           * point type.
           *
           * If timestep storage is fixed then time will have a scalar dataspace
           * and will contain the floating point increment between timesteps for
           * values stored.  It may have an attribute, offset, containing the
           * first floating point timestep value sampled.
           *
           * This dataset is optional.
           */
          H5::DataSet time;
      };

    }
  }
}

#endif  // HEMELB_IO_H5MD_TIMEDATA_H
