//
//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_IO_H5MD_FIXEDDATA_H
#define HEMELB_IO_H5MD_FIXEDATA_H

#include <cstdlib>
#include <string>
#include <atomic>

#include <boost/multi_array.hpp>

#include <hdf5.h>

#include "io/hdf5/H5Error.h"
#include "io/hdf5/H5Types.h"

namespace hemelb
{
  namespace io
  {
    namespace h5md
    {

      /**
       * Time-independent H5MD element.
       *
       * Contains a dataset "value" of type enumeration, string, float or
       * integer (H5::EnumType or subclasses of H5::AtomType).  Rank and
       * dimensionality of the data are fixed at creation time.
       */
      template <class ValueType, int Rank>
      class FixedData
      {
        public:

          static_assert(std::is_enum<ValueType>::value ||
                        std::is_integral<ValueType>::value ||
                        std::is_floating_point<ValueType>::value ||
                        std::is_assignable<ValueType, std::string>::value,
                        "ValueType must be one of enum, integer, floating point"
                        " or string");

          typedef ValueType value_type;

          typedef hsize_t size_type;

          static constexpr int rank = Rank;

          /**
           * Creates a H5MD Time-Dependent Data group with explicit storage of
           * timesteps.
           *
           * Timesteps are stored in a dataset separate from the samples.  A
           * "step" dataset of integer type is always created to store the
           * timestep indices.  An optional dataset "time" of floating point
           * type may be created to additionally store the timestep values.
           *
           * @param ExtentList must model the Collection concept
           * (http://www.boost.org/doc/libs/1_58_0/libs/utility/Collection.html)
           * ElementList::value_type must be assignable to hsize_t (unsigned
           * long long).
           *
           * @param location    where to create the group (typically a particles
           *                    or observables group)
           * @param name        the name of the group
           * @param dims        the dimensions of the data.  May have one or
           *                    more unlimited dimensions.
           * @param chunk_dims  if any dimension in dims is H5S_UNLIMITED then
           *                    the dataset must be stored in chunks.  This
           *                    specifies the size of each chunk.  See HDF5
           *                    docs for more details.
           * @param time        whether to create a "time" dataset to store
           *                    timestep values (default true)
           * @return a Time-Dependent Data group created at the specified
           *   location.
           */
          template <class ExtentList>
          static FixedData CreateExplicit(const hid_t & location,
                                         const std::string & name,
                                         const ExtentList & dims,
                                         const hid_t & hdf5_value_type = hdf5::Types<ValueType>::NATIVE) {
              assert(H5Iget_type(location) == H5I_GROUP);

              hid_t error;
              hid_t id;

              // Create a group with the specified name in the specified location
              HEMELB_HDF5_CALL(H5Gcreate, id, location, name.c_str(),
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
              std::shared_ptr<hid_t> group(new hid_t(id), hdf5::H5GroupDeleter);

              // Copy/convert the dimensions to an array of hsize_t expected by
              // HDF5
              hsize_t value_dims[rank];
              value_dims[0] = 0;
              std::copy(std::begin(dims), std::end(dims), std::begin(value_dims));

              // Create the value dataset
              HEMELB_HDF5_CALL(H5Screate_simple, id, rank,
                               value_dims, nullptr);
              std::shared_ptr<hid_t> value_space(new hid_t(id), hdf5::H5DataSpaceDeleter);

              HEMELB_HDF5_CALL(H5Dcreate, id, *group, "value",
                               hdf5_value_type, *value_space,
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
              std::shared_ptr<hid_t> value(new hid_t(id), hdf5::H5DataSetDeleter);

              return FixedData(group, name, value);
          }

          static FixedData Open(const hid_t & location, const std::string & name) {
            assert(H5Iget_type(location) == H5I_GROUP);

            hid_t id;

            // Open the group
            HEMELB_HDF5_CALL(H5Gopen, id, location, name.c_str(), H5P_DEFAULT);
            std::shared_ptr<hid_t> group(new hid_t(id), hdf5::H5GroupDeleter);

            // Open the value dataset and check the dataspace
            HEMELB_HDF5_CALL(H5Dopen, id, *group, "value", H5P_DEFAULT);
            std::shared_ptr<hid_t> value(new hid_t(id), hdf5::H5DataSetDeleter);
            HEMELB_HDF5_CALL(H5Dget_space, id, *value);
            std::shared_ptr<hid_t> value_space(new hid_t(id), hdf5::H5DataSpaceDeleter);
            assert(H5Sis_simple(*value_space) > 0);

            return FixedData(group, name, value);
          }

          template <class Sample, class = typename std::enable_if<Sample::dimensionality == Rank>::type>
          void SetData(const Sample & sample)
          {
          }

          template <class Sample, class IndexList,
                    class = typename std::enable_if<Sample::dimensionality == Rank &&
                                                    IndexList::dimensionality == Rank + 1>::type>
          void GetData(Sample & sample) const
          {
          }

        private:
          FixedData(const std::shared_ptr<hid_t> & group, const std::string & name,
                   const std::shared_ptr<hid_t> & value) :
              group(group), name(name), value(value)
          {
          }

          /** This group */
          std::shared_ptr<hid_t> group;

          /** The name of this group */
          std::string name;

          /**
           * Contains actual data values.
           */
          std::shared_ptr<hid_t> value;

      };

    }
  }
}

#endif  // HEMELB_IO_H5MD_FIXEDDATA_H
