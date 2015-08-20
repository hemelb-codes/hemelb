//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_IO_H5MD_TIMEDATA_H
#define HEMELB_IO_H5MD_TIMEDATA_H

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
       * Time-dependent H5MD element.
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
       * If time is explicitly stored then it will also contain a float-valued
       * dataset of the same rank and dimensionality as step, containing the
       * absolute values of the timesteps sampled.
       */
      template <class ValueType, int Rank,
                class StepType = unsigned long, class TimeType = double>
      class TimeData
      {
        public:

          static_assert(std::is_integral<StepType>::value,
                        "StepType must be integral");

          static_assert(std::is_floating_point<TimeType>::value,
                        "TimeType must be floating point");

          static_assert(std::is_enum<ValueType>::value ||
                        std::is_integral<ValueType>::value ||
                        std::is_floating_point<ValueType>::value ||
                        std::is_assignable<ValueType, std::string>::value,
                        "ValueType must be one of enum, integer, floating point"
                        " or string");

          typedef ValueType value_type;
          typedef StepType step_type;
          typedef TimeType time_type;

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
          static TimeData CreateExplicit(const hid_t & location,
                                         const std::string & name,
                                         const ExtentList & dims,
                                         bool time = true,
                                         const hid_t & hdf5_value_type = hdf5::Types<ValueType>::NATIVE,
                                         const hid_t & hdf5_step_type = hdf5::Types<StepType>::NATIVE,
                                         const hid_t & hdf5_time_type = hdf5::Types<TimeType>::NATIVE) {
              return CreateExplicit(location, name, dims, dims, time,
                                    hdf5_value_type, hdf5_step_type, hdf5_time_type);
          }
          template <class ExtentList>
          static TimeData CreateExplicit(const hid_t & location,
                                         const std::string & name,
                                         const ExtentList & dims,
                                         const ExtentList & chunk_dims,
                                         bool time = true,
                                         const hid_t & hdf5_value_type = hdf5::Types<ValueType>::NATIVE,
                                         const hid_t & hdf5_step_type = hdf5::Types<StepType>::NATIVE,
                                         const hid_t & hdf5_time_type = hdf5::Types<TimeType>::NATIVE) {
              assert(H5Iget_type(location) == H5I_GROUP);

              hid_t error;
              hid_t id;

              // Create a group with the specified name in the specified location
              HEMELB_HDF5_CALL(H5Gcreate, id, location, name.c_str(),
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
              std::shared_ptr<hid_t> group(new hid_t(id), hdf5::H5GroupDeleter);

              // Copy/convert the dimensions to an array of hsize_t expected by
              // HDF5 and prepend an H5S_UNLIMITED dimension for the timesteps
              hsize_t value_dims[rank + 1];
              value_dims[0] = 0;
              std::copy(std::begin(dims), std::end(dims), &value_dims[1]);
              hsize_t value_max_dims[rank + 1];
              value_max_dims[0] = H5S_UNLIMITED;
              std::copy(std::begin(dims), std::end(dims), &value_max_dims[1]);
              hsize_t value_chunk_dims[rank + 1];
              value_chunk_dims[0] = 1;
              std::copy(std::begin(chunk_dims), std::end(chunk_dims),
                        &value_chunk_dims[1]);

              // Create the value dataset
              HEMELB_HDF5_CALL(H5Screate_simple, id, rank + 1,
                               value_dims, value_max_dims);
              std::shared_ptr<hid_t> value_space(new hid_t(id), hdf5::H5DataSpaceDeleter);

              HEMELB_HDF5_CALL(H5Pcreate, id, H5P_DATASET_CREATE);
              std::shared_ptr<hid_t> create_plist(new hid_t(id), hdf5::H5PropertyListDeleter);
              HEMELB_HDF5_CALL(H5Pset_layout, error, *create_plist, H5D_CHUNKED);
              HEMELB_HDF5_CALL(H5Pset_chunk, error, *create_plist, rank + 1, value_chunk_dims);
              HEMELB_HDF5_CALL(H5Dcreate, id, *group, "value",
                               hdf5_value_type, *value_space,
                               H5P_DEFAULT, *create_plist, H5P_DEFAULT);
              std::shared_ptr<hid_t> value(new hid_t(id), hdf5::H5DataSetDeleter);

              // Create the timestep datasets
              const hsize_t timestep_dims[1] = { 0 };
              const hsize_t timestep_max_dims[1] = { H5S_UNLIMITED };
              const hsize_t timestep_chunk_dims[1] = { 1 };
              HEMELB_HDF5_CALL(H5Screate_simple, id, 1, timestep_dims, timestep_max_dims);
              std::shared_ptr<hid_t> timestep_space(new hid_t(id), hdf5::H5DataSpaceDeleter);

              HEMELB_HDF5_CALL(H5Pset_chunk, error, *create_plist, 1, timestep_chunk_dims);
              HEMELB_HDF5_CALL(H5Dcreate, id, *group, "step",
                               hdf5_step_type, *timestep_space,
                               H5P_DEFAULT, *create_plist, H5P_DEFAULT);
              std::shared_ptr<hid_t> step(new hid_t(id), hdf5::H5DataSetDeleter);

              if (time)
              {
                HEMELB_HDF5_CALL(H5Dcreate, id, *group, "time",
                                 hdf5_time_type, *timestep_space,
                                 H5P_DEFAULT, *create_plist, H5P_DEFAULT);
                std::shared_ptr<hid_t> time(new hid_t(id), hdf5::H5DataSetDeleter);
                return TimeData(group, name, value, step, time);
              }

              return TimeData(group, name, value, step);
          }

          /**
           * Creates a H5MD Time-Dependent Data group with timesteps stored
           * using a stride and offset.
           *
           * Timesteps are not stored explicitly in a dataset separate from the
           * samples.  Instead, the "step" and optional "time" datasets have
           * scalar dataspaces indicating the step and time increments,
           * respectively.  They may also have an "offset" attribute attached
           * indicating the first timestep stored.  If the floating point valued
           * timestep increment is zero or any of the offset attributes are zero
           * then they will not be created.
           *
           * @param ExtentList must model the Collection concept
           * (http://www.boost.org/doc/libs/1_58_0/libs/utility/Collection.html)
           * ElementList::value_type must be assignable to hsize_t (unsigned
           * long long).
           *
           * @param location        where to create the group (typically a particles
           *                        or observables group)
           * @param name            the name of the group
           * @param dims            the dimensions of the data.  May have one or
           *                        more unlimited dimensions.
           * @param chunk_dims      if any dimension in dims is H5S_UNLIMITED then
           *                        the dataset must be stored in chunks.  This
           *                        specifies the size of each chunk.  See HDF5
           *                        docs for more details.
           * @param time            whether to create a "time" dataset to store
           *                        timestep values (default true)
           * @param step_increment  the integer increment between timesteps
           *                        (default 1)
           * @param step_offset     the first integer timestep (default 0)
           * @param time_increment  the floating point increment between
           *                        timesteps (default 0.0)
           * @param time_offset     the first floating point timestep (default
           *                        0.0)
           * @return a Time-Dependent Data group created at the specified
           *   location.
           */
          template <class ExtentList>
          static TimeData CreateFixed(const hid_t & location,
                                      const std::string & name,
                                      const ExtentList & dims,
                                      bool time = true,
                                      const StepType & step_increment = 1,
                                      const StepType & step_offset = 0,
                                      const TimeType & time_increment = 0.0,
                                      const TimeType & time_offset = 0.0,
                                      const hid_t & hdf5_value_type = hdf5::Types<ValueType>::NATIVE,
                                      const hid_t & hdf5_step_type = hdf5::Types<StepType>::NATIVE,
                                      const hid_t & hdf5_time_type = hdf5::Types<TimeType>::NATIVE) {
              return CreateFixed(location, name, dims, dims, time, step_increment,
                                 step_offset, time_increment, time_offset,
                                 hdf5_value_type, hdf5_step_type, hdf5_time_type);
          }
          template <class ExtentList>
          static TimeData CreateFixed(const hid_t & location,
                                      const std::string & name,
                                      const ExtentList & dims,
                                      const ExtentList & chunk_dims,
                                      bool time = true,
                                      const StepType & step_increment = 1,
                                      const StepType & step_offset = 0,
                                      const TimeType & time_increment = 0.0,
                                      const TimeType & time_offset = 0.0,
                                      const hid_t & hdf5_value_type = hdf5::Types<ValueType>::NATIVE,
                                      const hid_t & hdf5_step_type = hdf5::Types<StepType>::NATIVE,
                                      const hid_t & hdf5_time_type = hdf5::Types<TimeType>::NATIVE) {
              assert(H5Iget_type(location) == H5I_GROUP);

              hid_t error;
              hid_t id;

              // Create a group with the specified name in the specified location
              HEMELB_HDF5_CALL(H5Gcreate, id, location, name.c_str(),
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
              std::shared_ptr<hid_t> group(new hid_t(id), hdf5::H5GroupDeleter);

              // Copy/convert the dimensions to an array of hsize_t expected by
              // HDF5 and prepend an H5S_UNLIMITED dimension for the timesteps
              hsize_t value_dims[rank + 1];
              value_dims[0] = 0;
              std::copy(std::begin(dims), std::end(dims), &value_dims[1]);
              hsize_t value_max_dims[rank + 1];
              value_max_dims[0] = H5S_UNLIMITED;
              std::copy(std::begin(dims), std::end(dims), &value_max_dims[1]);
              hsize_t value_chunk_dims[rank + 1];
              value_chunk_dims[0] = 1;
              std::copy(std::begin(chunk_dims), std::end(chunk_dims),
                        &value_chunk_dims[1]);

              // Create the value dataset
              HEMELB_HDF5_CALL(H5Screate_simple, id, rank + 1,
                               value_dims, value_max_dims);
              std::shared_ptr<hid_t> value_space(new hid_t(id), hdf5::H5DataSpaceDeleter);

              HEMELB_HDF5_CALL(H5Pcreate, id, H5P_DATASET_CREATE);
              std::shared_ptr<hid_t> create_plist(new hid_t(id), hdf5::H5PropertyListDeleter);
              HEMELB_HDF5_CALL(H5Pset_layout, error, *create_plist, H5D_CHUNKED);
              HEMELB_HDF5_CALL(H5Pset_chunk, error, *create_plist, rank + 1, value_chunk_dims);
              HEMELB_HDF5_CALL(H5Dcreate, id, *group, "value",
                               hdf5_value_type, *value_space,
                               H5P_DEFAULT, *create_plist, H5P_DEFAULT);
              std::shared_ptr<hid_t> value(new hid_t(id), hdf5::H5DataSetDeleter);

              // Create the scalar timestep datasets
              HEMELB_HDF5_CALL(H5Screate, id, H5S_SCALAR);
              std::shared_ptr<hid_t> timestep_space(new hid_t(id), hdf5::H5DataSpaceDeleter);
              HEMELB_HDF5_CALL(H5Dcreate, id, *group, "step",
                               hdf5_step_type, *timestep_space,
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
              std::shared_ptr<hid_t> step(new hid_t(id), hdf5::H5DataSetDeleter);

              // Attach the offset attribute
              HEMELB_HDF5_CALL(H5Acreate, id, *step, "offset",
                               hdf5_step_type, *timestep_space,
                               H5P_DEFAULT, H5P_DEFAULT);
              std::shared_ptr<hid_t> step_offset_attr(new hid_t(id), hdf5::H5AttributeDeleter);
              HEMELB_HDF5_CALL(H5Awrite, error, *step_offset_attr,
                               hdf5::Types<step_type>::NATIVE, &step_offset);

              if (time)
              {
                HEMELB_HDF5_CALL(H5Dcreate, id, *group, "time",
                                 hdf5_time_type, *timestep_space,
                                 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                std::shared_ptr<hid_t> time(new hid_t(id), hdf5::H5DataSetDeleter);

                HEMELB_HDF5_CALL(H5Acreate, id, *time, "offset",
                                 hdf5_time_type, *timestep_space,
                                 H5P_DEFAULT, H5P_DEFAULT);
                std::shared_ptr<hid_t> time_offset_attr(new hid_t(id), hdf5::H5AttributeDeleter);

                HEMELB_HDF5_CALL(H5Awrite, error, *time_offset_attr,
                                 hdf5_time_type, &time_offset);

                return TimeData(group, name, value, step, time);
              }

              return TimeData(group, name, value, step);
          }

          static TimeData Open(const hid_t & location, const std::string & name) {
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

            // Open the timestep datasets
            HEMELB_HDF5_CALL(H5Dopen, id, *group, "step", H5P_DEFAULT);
            std::shared_ptr<hid_t> step(new hid_t(id), hdf5::H5DataSetDeleter);

            std::shared_ptr<hid_t> time;
            if ((id = H5Dopen(*group, "time", H5P_DEFAULT)) >= 0)
            {
              time.reset(new hid_t(id), hdf5::H5DataSetDeleter);
            }

            return TimeData(group, name, value, step, time);
          }

          bool IsFixed() const {
            // Fixed timestep data has a scalar dataspace for the step and time
            // datasets.
            // Explicit timestep data has a simple dataspace for the step and
            // time datasets.
            // Only the step dataset is guaranteed to exist - time may not (and
            // value is not useful here).
            hid_t id;

            HEMELB_HDF5_CALL(H5Dget_space, id, *step);
            std::shared_ptr<hid_t> step_space(new hid_t(id), hdf5::H5DataSpaceDeleter);

            // get the rank of the dataspace
            int rank;
            HEMELB_HDF5_CALL(H5Sget_simple_extent_ndims, rank, *step_space);

            // scalar dataspaces have rank 0
            return rank == 0;
          }

          bool HasTime() const
          {
            return time;
          }

          size_type GetNumSamples() const
          {
            herr_t error;
            hid_t id;

            HEMELB_HDF5_CALL(H5Dget_space, id, *value);
            std::shared_ptr<hid_t> space(new hid_t(id), hdf5::H5DataSpaceDeleter);

            hsize_t dims[rank + 1];
            HEMELB_HDF5_CALL(H5Sget_simple_extent_dims, error, *space, dims, nullptr);

            return dims[0];
          }

          template <class Sample, class = typename std::enable_if<Sample::dimensionality == Rank>::type>
          void AddSample(const Sample & sample, const step_type & step, const time_type & time)
          {
              herr_t error;
              hid_t id;

              // Get the current extents of the dataset
              HEMELB_HDF5_CALL(H5Dget_space, id, *value);
              std::shared_ptr<hid_t> value_file_space(new hid_t(id), hdf5::H5DataSpaceDeleter);
              hsize_t dims[rank + 1], max_dims[rank + 1];
              HEMELB_HDF5_CALL(H5Sget_simple_extent_dims, error, *value_file_space, dims, max_dims);

              // Expand the extents of the dataset
              ++dims[0];
              HEMELB_HDF5_CALL(H5Dset_extent, error, *value, dims);

              // Refresh the dataspace
              HEMELB_HDF5_CALL(H5Dget_space, id, *value);
              value_file_space.reset(new hid_t(id), hdf5::H5DataSpaceDeleter);

              // Create a dataspace to represent the extension to the dataset
              hsize_t ext[rank + 1];
              ext[0] = 1;
              std::copy(&dims[1], &dims[rank + 1], &ext[1]);
              HEMELB_HDF5_CALL(H5Screate_simple, id, rank, ext, nullptr);
              std::shared_ptr<hid_t> value_mem_space(new hid_t(id), hdf5::H5DataSpaceDeleter);

              // Select a hyperslab in the extended part of the dataset
              --dims[0];
              HEMELB_HDF5_CALL(H5Sselect_hyperslab, error, *value_file_space,
                               H5S_SELECT_SET, dims, nullptr, ext, nullptr);

              // Write the data to the extended part of the dataset
              HEMELB_HDF5_CALL(H5Dwrite, error, *value, hdf5::Types<value_type>::NATIVE,
                               *value_mem_space, *value_file_space, H5P_DEFAULT, sample.data());
          }

          template <class Sample, class IndexList,
                    class = typename std::enable_if<Sample::dimensionality == Rank &&
                                                    IndexList::dimensionality == Rank + 1>::type>
          void GetSample(const IndexList & index, Sample & sample) const
          {
              herr_t error;
              hid_t id;

              // Get the dataspace from the file
              HEMELB_HDF5_CALL(H5Dget_space, id, *value);
              std::shared_ptr<hid_t> file_space(new hid_t(id), hdf5::H5DataSpaceDeleter);

              // Select a hyperslab in the dataset representing the sample in
              // memory
              const hsize_t slab[IndexList::dimensionality];
              std::copy(std::begin(index), std::end(index), std::begin(slab));
              const hsize_t count[IndexList::dimensionality];
              std::fill(std::begin(count), std::end(count), 1);
              HEMELB_HDF5_CALL(H5Sselect_hyperslab, error, *file_space,
                               H5S_SELECT_SET, slab, nullptr, count, nullptr);

              // Create a dataspace representing the value in memory
              hsize_t dims[Sample::dimensionality];
              std::copy(std::begin(index) + 1, std::end(index), std::begin(dims));
              HEMELB_HDF5_CALL(H5Screate_simple, id, Sample::dimensionality, dims, nullptr);
              std::shared_ptr<hid_t> mem_space(new hid_t(id), hdf5::H5DataSpaceDeleter);

              // Read the data from the slab
              HEMELB_HDF5_CALL(H5Dread, error, *value, hdf5::Types<typename Sample::value_type>::NATIVE,
                               *mem_space, *file_space, H5P_DEFAULT, sample.data());
          }

          step_type GetStep(const size_type & index) const
          {
            step_type res = 0;

            herr_t error;
            hid_t id;

            // Get the dataspace from the file
            HEMELB_HDF5_CALL(H5Dget_space, id, *step);
            std::shared_ptr<hid_t> file_space(new hid_t(id), hdf5::H5DataSpaceDeleter);

            // Create a scalar dataspace representing the value in memory
            HEMELB_HDF5_CALL(H5Screate, id, H5S_SCALAR);
            std::shared_ptr<hid_t> mem_space(new hid_t(id), hdf5::H5DataSpaceDeleter);

            // If the dataspace is simple then it contains all the timesteps
            // in ascending order
            htri_t simple;
            HEMELB_HDF5_CALL(H5Sis_simple, simple, *file_space);
            if (simple > 0)
            {
              // Select a hyperslab in the dataset
              const hsize_t slab[1] = { index };
              const hsize_t count[1] = { 1 };
              HEMELB_HDF5_CALL(H5Sselect_hyperslab, error, *file_space,
                               H5S_SELECT_SET, slab, nullptr, count, nullptr);

              // Read the data from the hyperslab
              HEMELB_HDF5_CALL(H5Dread, error, *step, hdf5::Types<step_type>::NATIVE,
                               *mem_space, *file_space, H5P_DEFAULT, &res);
            }
            else  // The dataspace is scalar and contains the increment
            {
              step_type increment, offset = 0;

              // Read the increment from the dataset
              HEMELB_HDF5_CALL(H5Dread, error, *step, hdf5::Types<step_type>::NATIVE,
                               *mem_space, *file_space, H5P_DEFAULT, &increment);

              // If there is an offset attribute read that too
              htri_t has_offset_attr;
              HEMELB_HDF5_CALL(H5Aexists, has_offset_attr, *step, "offset");
              if (has_offset_attr > 0)
              {
                // Open the offset attribute
                HEMELB_HDF5_CALL(H5Aopen, id, *step, "offset", H5P_DEFAULT);
                std::shared_ptr<hid_t> offset_attr(new hid_t(id), hdf5::H5AttributeDeleter);

                // Read the offset
                HEMELB_HDF5_CALL(H5Aread, error, *offset_attr,
                                 hdf5::Types<step_type>::NATIVE, &offset);
              }

              // Multiply by the index and add the offset
              res = static_cast<step_type>(index) * increment + offset;
            }

            return res;
          }

          time_type GetTime(const size_type & index) const
          {
            time_type res = 0.0;

            if (time)
            {
              herr_t error;
              hid_t id;

              // Get the dataspace from the file
              HEMELB_HDF5_CALL(H5Dget_space, id, *time);
              std::shared_ptr<hid_t> file_space(new hid_t(id), hdf5::H5DataSpaceDeleter);

              // Create a scalar dataspace representing the value in memory
              HEMELB_HDF5_CALL(H5Screate, id, H5S_SCALAR);
              std::shared_ptr<hid_t> mem_space(new hid_t(id), hdf5::H5DataSpaceDeleter);

              // If the dataspace is simple then it contains all the timesteps
              // in ascending order
              htri_t simple;
              HEMELB_HDF5_CALL(H5Sis_simple, simple, *file_space);
              if (simple > 0)
              {
                // Select a hyperslab in the dataset
                const hsize_t slab[1] = { index };
                const hsize_t count[1] = { 1 };
                HEMELB_HDF5_CALL(H5Sselect_hyperslab, error, *file_space,
                                 H5S_SELECT_SET, slab, nullptr, count, nullptr);

                // Read the data from the hyperslab
                HEMELB_HDF5_CALL(H5Dread, error, *time, hdf5::Types<time_type>::NATIVE,
                                 *mem_space, *file_space, H5P_DEFAULT, &res);
              }
              else  // The dataspace is scalar and contains the increment
              {
                time_type increment, offset = 0.0;

                // Read the increment from the dataset
                HEMELB_HDF5_CALL(H5Dread, error, *time, hdf5::Types<time_type>::NATIVE,
                                 *mem_space, *file_space, H5P_DEFAULT, &increment);

                // If there is an offset attribute read that too
                htri_t has_offset_attr;
                HEMELB_HDF5_CALL(H5Aexists, has_offset_attr, *time, "offset");
                if (has_offset_attr > 0)
                {
                  // Open the offset attribute
                  HEMELB_HDF5_CALL(H5Aopen, id, *time, "offset", H5P_DEFAULT);
                  std::shared_ptr<hid_t> offset_attr(new hid_t(id), hdf5::H5AttributeDeleter);

                  // Read the offset
                  HEMELB_HDF5_CALL(H5Aread, error, *offset_attr,
                                   hdf5::Types<time_type>::NATIVE, &offset);
                }

                // Multiply by the index and add the offset
                res = static_cast<time_type>(index) * increment + offset;
              }
            }

            return res;
          }

        private:
          TimeData(const std::shared_ptr<hid_t> & group, const std::string & name,
                   const std::shared_ptr<hid_t> & value,
                   const std::shared_ptr<hid_t> & step,
                   const std::shared_ptr<hid_t> & time) :
              group(group), name(name), value(value), step(step), time(time)
          {
          }
          TimeData(const std::shared_ptr<hid_t> & group, const std::string & name,
                   const std::shared_ptr<hid_t> & value,
                   const std::shared_ptr<hid_t> & step) :
              group(group), name(name), value(value), step(step)
          {
          }

          /** This group */
          std::shared_ptr<hid_t> group;

          /** The name of this group */
          std::string name;

          /**
           * Contains actual data values from the simulation.
           *
           * Rank is one plus the dimensionality of the data to account for the
           * timestep dimension.  The timestep dimension is unlimited.
           */
          std::shared_ptr<hid_t> value;

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
          std::shared_ptr<hid_t> step;

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
          std::shared_ptr<hid_t> time;
      };

    }
  }
}

#endif  // HEMELB_IO_H5MD_TIMEDATA_H
