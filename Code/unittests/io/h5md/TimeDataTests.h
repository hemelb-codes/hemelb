//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#include "io/h5md/TimeData.h"

#include <cppunit/extensions/HelperMacros.h>
#include <hdf5.h>

#include "unittests/helpers/FolderTestFixture.h"

#include "io/hdf5/H5Error.h"

#ifndef HEMELB_UNITTESTS_IO_H5MD_TIMEDATATESTS_H
#define HEMELB_UNITTESTS_IO_H5MD_TIMEDATATESTS_H

namespace hemelb
{
  namespace unittests
  {
    namespace io
    {
      namespace h5md
      {

        class TimeDataTests : public hemelb::unittests::helpers::FolderTestFixture
        {
            CPPUNIT_TEST_SUITE(TimeDataTests);
            CPPUNIT_TEST(testCreateExplicit);
            CPPUNIT_TEST(testCreateExplicitNoTime);
            CPPUNIT_TEST(testCreateFixed);
            CPPUNIT_TEST(testCreateFixedNoTime);
            CPPUNIT_TEST_SUITE_END();

          public:
            void setUp()
            {
              herr_t error;

              // Disable HDF5's automatic printing of error messages
              HEMELB_HDF5_CALL(H5Eset_auto2, error, H5E_DEFAULT, nullptr, nullptr);

              // Use the in-memory "core" driver for unit testing
              hid_t fapl;
              HEMELB_HDF5_CALL(H5Pcreate, fapl, H5P_FILE_ACCESS);
              HEMELB_HDF5_CALL(H5Pset_fapl_core, error, fapl, 1024 * 1024, false);

              // Create an empty HDF5 file
              HEMELB_HDF5_CALL(H5Fcreate, file, "file.h5md", H5F_ACC_EXCL, H5P_DEFAULT, fapl);

              // Create a "particles" group for testing
              HEMELB_HDF5_CALL(H5Gcreate, location, file, "particles", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

              // Close the file access property list
              HEMELB_HDF5_CALL(H5Pclose, error, fapl);
            }

            void tearDown()
            {
              herr_t error;
              HEMELB_HDF5_CALL(H5Gclose, error, location);
              HEMELB_HDF5_CALL(H5Fclose, error, file);
            }

            void testCreateExplicit() {
              // Create a rank-1, 3-dimensional, time-dependent dataset with
              // explicit step and time storage
              typedef hemelb::io::h5md::TimeData<double, 1> TimeData;
              const hsize_t dims[] = { 3 };
              TimeData timedata = TimeData::CreateExplicit(location, "RBC", dims);

              CPPUNIT_ASSERT(!timedata.IsFixed());

              hid_t id;
              herr_t error;

              // A group exists with the correct name at the correct location
              id = H5Gopen(location, "RBC", H5P_DEFAULT);
              CPPUNIT_ASSERT_MESSAGE("time dependent group doesn't exist", id >= 0);
              std::shared_ptr<hid_t> group(new hid_t(id), hemelb::io::hdf5::H5GroupDeleter);

              // The value dataset exists within the group
              id = H5Dopen(*group, "value", H5P_DEFAULT);
              CPPUNIT_ASSERT_MESSAGE("value dataset doesn't exist", id >= 0);
              std::shared_ptr<hid_t> value(new hid_t(id), hemelb::io::hdf5::H5DataSetDeleter);

              HEMELB_HDF5_CALL(H5Dget_space, id, *value);
              std::shared_ptr<hid_t> value_space(new hid_t(id), hemelb::io::hdf5::H5DataSpaceDeleter);
              CPPUNIT_ASSERT_MESSAGE("value dataspace is not simple",
                                     H5Sis_simple(*value_space) > 0);
              const int value_rank = H5Sget_simple_extent_ndims(*value_space);
              CPPUNIT_ASSERT_MESSAGE("value dataspace rank is not 1 plus the "
                                     "tensor rank of the data stored",
                                     value_rank == 2);
              hsize_t value_shape[value_rank];
              HEMELB_HDF5_CALL(H5Sget_simple_extent_dims, error, *value_space,
                               nullptr, value_shape);
              CPPUNIT_ASSERT_MESSAGE("value shape is not the shape of a single data item",
                                     value_shape[1] == 3);
              CPPUNIT_ASSERT_MESSAGE("value shape is not prepended by a [variable] dimension",
                                     value_shape[0] == H5S_UNLIMITED);

              // The step dataset exists within the group
              id = H5Dopen(*group, "step", H5P_DEFAULT);
              CPPUNIT_ASSERT_MESSAGE("step dataset doesn't exist", id >= 0);
              std::shared_ptr<hid_t> step(new hid_t(id), hemelb::io::hdf5::H5DataSetDeleter);

              CPPUNIT_ASSERT_MESSAGE("step type is not integer",
                                     H5Tget_class(H5Dget_type(*step)) == H5T_INTEGER);

              HEMELB_HDF5_CALL(H5Dget_space, id, *step);
              std::shared_ptr<hid_t> step_space(new hid_t(id), hemelb::io::hdf5::H5DataSpaceDeleter);
              CPPUNIT_ASSERT_MESSAGE("step dataspace is not simple",
                                     H5Sis_simple(*step_space) > 0);
              const int step_rank = H5Sget_simple_extent_ndims(*step_space);
              CPPUNIT_ASSERT_MESSAGE("step dataspace rank is not 1", step_rank == 1);
              hsize_t step_shape[step_rank];
              HEMELB_HDF5_CALL(H5Sget_simple_extent_dims, error, *step_space,
                               nullptr, step_shape);
              CPPUNIT_ASSERT_MESSAGE("step shape is not [variable]",
                                     step_shape[0] == H5S_UNLIMITED);

              // The time dataset exists within the group
              id = H5Dopen(*group, "time", H5P_DEFAULT);
              CPPUNIT_ASSERT_MESSAGE("time dataset doesn't exist", id >= 0);
              std::shared_ptr<hid_t> time(new hid_t(id), hemelb::io::hdf5::H5DataSetDeleter);

              CPPUNIT_ASSERT_MESSAGE("time type is not floating point or integer",
                                     H5Tget_class(H5Dget_type(*time)) == H5T_FLOAT ||
                                     H5Tget_class(H5Dget_type(*time)) == H5T_INTEGER);

              HEMELB_HDF5_CALL(H5Dget_space, id, *time);
              std::shared_ptr<hid_t> time_space(new hid_t(id), hemelb::io::hdf5::H5DataSpaceDeleter);
              CPPUNIT_ASSERT_MESSAGE("time dataspace is not simple",
                                     H5Sis_simple(*time_space) > 0);
              const int time_rank = H5Sget_simple_extent_ndims(*time_space);
              CPPUNIT_ASSERT_MESSAGE("time dataspace rank is not 1", time_rank == 1);
              hsize_t time_shape[time_rank];
              HEMELB_HDF5_CALL(H5Sget_simple_extent_dims, error, *time_space,
                               nullptr, time_shape);
              CPPUNIT_ASSERT_MESSAGE("time shape is not [variable]",
                                     time_shape[0] == H5S_UNLIMITED);
            }

            void testCreateExplicitNoTime() {
              // Create a rank-1, 3-dimensional, time-dependent dataset with
              // explicit step and time storage
              typedef hemelb::io::h5md::TimeData<double, 1> TimeData;
              const hsize_t dims[] = { 3 };
              TimeData timedata = TimeData::CreateExplicit(location, "RBC", dims, dims, false);

              CPPUNIT_ASSERT(!timedata.IsFixed());

              hid_t id;
              herr_t error;

              // A group exists with the correct name at the correct location
              id = H5Gopen(location, "RBC", H5P_DEFAULT);
              CPPUNIT_ASSERT_MESSAGE("time dependent group doesn't exist", id >= 0);
              std::shared_ptr<hid_t> group(new hid_t(id), hemelb::io::hdf5::H5GroupDeleter);

              // The value dataset exists within the group
              id = H5Dopen(*group, "value", H5P_DEFAULT);
              CPPUNIT_ASSERT_MESSAGE("value dataset doesn't exist", id >= 0);
              std::shared_ptr<hid_t> value(new hid_t(id), hemelb::io::hdf5::H5DataSetDeleter);

              HEMELB_HDF5_CALL(H5Dget_space, id, *value);
              std::shared_ptr<hid_t> value_space(new hid_t(id), hemelb::io::hdf5::H5DataSpaceDeleter);
              CPPUNIT_ASSERT_MESSAGE("value dataspace is not simple",
                                     H5Sis_simple(*value_space) > 0);
              const int value_rank = H5Sget_simple_extent_ndims(*value_space);
              CPPUNIT_ASSERT_MESSAGE("value dataspace rank is not 1 plus the "
                                     "tensor rank of the data stored",
                                     value_rank == 2);
              hsize_t value_shape[value_rank];
              HEMELB_HDF5_CALL(H5Sget_simple_extent_dims, error, *value_space,
                               nullptr, value_shape);
              CPPUNIT_ASSERT_MESSAGE("value shape is not the shape of a single data item",
                                     value_shape[1] == 3);
              CPPUNIT_ASSERT_MESSAGE("value shape is not prepended by a [variable] dimension",
                                     value_shape[0] == H5S_UNLIMITED);

              // The step dataset exists within the group
              id = H5Dopen(*group, "step", H5P_DEFAULT);
              CPPUNIT_ASSERT_MESSAGE("step dataset doesn't exist", id >= 0);
              std::shared_ptr<hid_t> step(new hid_t(id), hemelb::io::hdf5::H5DataSetDeleter);

              CPPUNIT_ASSERT_MESSAGE("step type is not integer",
                                     H5Tget_class(H5Dget_type(*step)) == H5T_INTEGER);

              HEMELB_HDF5_CALL(H5Dget_space, id, *step);
              std::shared_ptr<hid_t> step_space(new hid_t(id), hemelb::io::hdf5::H5DataSpaceDeleter);
              CPPUNIT_ASSERT_MESSAGE("step dataspace is not simple",
                                     H5Sis_simple(*step_space) > 0);
              const int step_rank = H5Sget_simple_extent_ndims(*step_space);
              CPPUNIT_ASSERT_MESSAGE("step dataspace rank is not 1", step_rank == 1);
              hsize_t step_shape[step_rank];
              HEMELB_HDF5_CALL(H5Sget_simple_extent_dims, error, *step_space,
                               nullptr, step_shape);
              CPPUNIT_ASSERT_MESSAGE("step shape is not [variable]",
                                     step_shape[0] == H5S_UNLIMITED);

              // The time dataset doesn't exist within the group
              id = H5Dopen(*group, "time", H5P_DEFAULT);
              CPPUNIT_ASSERT_MESSAGE("time dataset exists", id < 0);
              if (id >= 0)
                H5Dclose(id);
            }

            void testCreateFixed() {
              // Create a rank-1, 3-dimensional, time-dependent dataset with
              // explicit step and time storage
              typedef hemelb::io::h5md::TimeData<double, 1> TimeData;
              const hsize_t dims[] = { 3 };
              TimeData timedata = TimeData::CreateFixed(location, "RBC", dims);

              CPPUNIT_ASSERT(timedata.IsFixed());

              hid_t id;
              herr_t error;

              // A group exists with the correct name at the correct location
              id = H5Gopen(location, "RBC", H5P_DEFAULT);
              CPPUNIT_ASSERT_MESSAGE("time dependent group doesn't exist", id >= 0);
              std::shared_ptr<hid_t> group(new hid_t(id), hemelb::io::hdf5::H5GroupDeleter);

              // The value dataset exists within the group
              id = H5Dopen(*group, "value", H5P_DEFAULT);
              CPPUNIT_ASSERT_MESSAGE("value dataset doesn't exist", id >= 0);
              std::shared_ptr<hid_t> value(new hid_t(id), hemelb::io::hdf5::H5DataSetDeleter);

              HEMELB_HDF5_CALL(H5Dget_space, id, *value);
              std::shared_ptr<hid_t> value_space(new hid_t(id), hemelb::io::hdf5::H5DataSpaceDeleter);
              CPPUNIT_ASSERT_MESSAGE("value dataspace is not simple",
                                     H5Sis_simple(*value_space) > 0);
              const int value_rank = H5Sget_simple_extent_ndims(*value_space);
              CPPUNIT_ASSERT_MESSAGE("value dataspace rank is not 1 plus the "
                                     "tensor rank of the data stored",
                                     value_rank == 2);
              hsize_t value_shape[value_rank];
              HEMELB_HDF5_CALL(H5Sget_simple_extent_dims, error, *value_space,
                               nullptr, value_shape);
              CPPUNIT_ASSERT_MESSAGE("value shape is not the shape of a single data item",
                                     value_shape[1] == 3);
              CPPUNIT_ASSERT_MESSAGE("value shape is not prepended by a [variable] dimension",
                                     value_shape[0] == H5S_UNLIMITED);

              // The step dataset exists within the group
              id = H5Dopen(*group, "step", H5P_DEFAULT);
              CPPUNIT_ASSERT_MESSAGE("step dataset doesn't exist", id >= 0);
              std::shared_ptr<hid_t> step(new hid_t(id), hemelb::io::hdf5::H5DataSetDeleter);

              CPPUNIT_ASSERT_MESSAGE("step type is not integer",
                                     H5Tget_class(H5Dget_type(*step)) == H5T_INTEGER);

              HEMELB_HDF5_CALL(H5Dget_space, id, *step);
              std::shared_ptr<hid_t> step_space(new hid_t(id), hemelb::io::hdf5::H5DataSpaceDeleter);
              int rank;
              HEMELB_HDF5_CALL(H5Sget_simple_extent_ndims, rank, *step_space);
              CPPUNIT_ASSERT_MESSAGE("step dataspace is not scalar", rank == 0);

              HEMELB_HDF5_CALL(H5Aopen, id, *step, "offset", H5P_DEFAULT);
              std::shared_ptr<hid_t> attr(new hid_t(id), hemelb::io::hdf5::H5AttributeDeleter);
              CPPUNIT_ASSERT_MESSAGE("step offset type is not integer",
                                     H5Tget_class(H5Aget_type(*attr)) == H5T_INTEGER);
              HEMELB_HDF5_CALL(H5Dget_space, id, *step);
              std::shared_ptr<hid_t> attr_space(new hid_t(id), hemelb::io::hdf5::H5DataSpaceDeleter);
              HEMELB_HDF5_CALL(H5Sget_simple_extent_ndims, rank, *attr_space);
              CPPUNIT_ASSERT_MESSAGE("step offset dataspace is not scalar",
                                     rank == 0);


              // The time dataset exists within the group
              id = H5Dopen(*group, "time", H5P_DEFAULT);
              CPPUNIT_ASSERT_MESSAGE("time dataset doesn't exist", id >= 0);
              std::shared_ptr<hid_t> time(new hid_t(id), hemelb::io::hdf5::H5DataSetDeleter);

              CPPUNIT_ASSERT_MESSAGE("time type is not floating point or integer",
                                     H5Tget_class(H5Dget_type(*time)) == H5T_FLOAT ||
                                     H5Tget_class(H5Dget_type(*time)) == H5T_INTEGER);

              HEMELB_HDF5_CALL(H5Dget_space, id, *time);
              std::shared_ptr<hid_t> time_space(new hid_t(id), hemelb::io::hdf5::H5DataSpaceDeleter);
              HEMELB_HDF5_CALL(H5Sget_simple_extent_ndims, rank, *time_space);
              CPPUNIT_ASSERT_MESSAGE("time dataspace is not scalar",
                                     rank == 0);

              HEMELB_HDF5_CALL(H5Aopen, id, *time, "offset", H5P_DEFAULT);
              attr.reset(new hid_t(id), hemelb::io::hdf5::H5AttributeDeleter);
              CPPUNIT_ASSERT_MESSAGE("time offset type is not floating point or integer",
                                     H5Tget_class(H5Aget_type(*attr)) == H5T_FLOAT ||
                                     H5Tget_class(H5Aget_type(*attr)) == H5T_INTEGER);
              HEMELB_HDF5_CALL(H5Dget_space, id, *step);
              attr_space.reset(new hid_t(id), hemelb::io::hdf5::H5DataSpaceDeleter);
              HEMELB_HDF5_CALL(H5Sget_simple_extent_ndims, rank, *attr_space);
              CPPUNIT_ASSERT_MESSAGE("time offset dataspace is not scalar",
                                     rank == 0);
            }

            void testCreateFixedNoTime() {
              // Create a rank-1, 3-dimensional, time-dependent dataset with
              // explicit step and time storage
              typedef hemelb::io::h5md::TimeData<double, 1> TimeData;
              const hsize_t dims[] = { 3 };
              TimeData timedata = TimeData::CreateFixed(location, "RBC", dims,
                                                        dims, false);

              CPPUNIT_ASSERT(timedata.IsFixed());

              hid_t id;
              herr_t error;

              // A group exists with the correct name at the correct location
              id = H5Gopen(location, "RBC", H5P_DEFAULT);
              CPPUNIT_ASSERT_MESSAGE("time dependent group doesn't exist", id >= 0);
              std::shared_ptr<hid_t> group(new hid_t(id), hemelb::io::hdf5::H5GroupDeleter);

              // The value dataset exists within the group
              id = H5Dopen(*group, "value", H5P_DEFAULT);
              CPPUNIT_ASSERT_MESSAGE("value dataset doesn't exist", id >= 0);
              std::shared_ptr<hid_t> value(new hid_t(id), hemelb::io::hdf5::H5DataSetDeleter);

              HEMELB_HDF5_CALL(H5Dget_space, id, *value);
              std::shared_ptr<hid_t> value_space(new hid_t(id), hemelb::io::hdf5::H5DataSpaceDeleter);
              CPPUNIT_ASSERT_MESSAGE("value dataspace is not simple",
                                     H5Sis_simple(*value_space) > 0);
              const int value_rank = H5Sget_simple_extent_ndims(*value_space);
              CPPUNIT_ASSERT_MESSAGE("value dataspace rank is not 1 plus the "
                                     "tensor rank of the data stored",
                                     value_rank == 2);
              hsize_t value_shape[value_rank];
              HEMELB_HDF5_CALL(H5Sget_simple_extent_dims, error, *value_space,
                               nullptr, value_shape);
              CPPUNIT_ASSERT_MESSAGE("value shape is not the shape of a single data item",
                                     value_shape[1] == 3);
              CPPUNIT_ASSERT_MESSAGE("value shape is not prepended by a [variable] dimension",
                                     value_shape[0] == H5S_UNLIMITED);

              // The step dataset exists within the group
              id = H5Dopen(*group, "step", H5P_DEFAULT);
              CPPUNIT_ASSERT_MESSAGE("step dataset doesn't exist", id >= 0);
              std::shared_ptr<hid_t> step(new hid_t(id), hemelb::io::hdf5::H5DataSetDeleter);

              CPPUNIT_ASSERT_MESSAGE("step type is not integer",
                                     H5Tget_class(H5Dget_type(*step)) == H5T_INTEGER);

              HEMELB_HDF5_CALL(H5Dget_space, id, *step);
              std::shared_ptr<hid_t> step_space(new hid_t(id), hemelb::io::hdf5::H5DataSpaceDeleter);
              int rank;
              HEMELB_HDF5_CALL(H5Sget_simple_extent_ndims, rank, *step_space);
              CPPUNIT_ASSERT_MESSAGE("step dataspace is not scalar", rank == 0);

              HEMELB_HDF5_CALL(H5Aopen, id, *step, "offset", H5P_DEFAULT);
              std::shared_ptr<hid_t> attr(new hid_t(id), hemelb::io::hdf5::H5AttributeDeleter);
              CPPUNIT_ASSERT_MESSAGE("step offset type is not integer",
                                     H5Tget_class(H5Aget_type(*attr)) == H5T_INTEGER);
              HEMELB_HDF5_CALL(H5Dget_space, id, *step);
              std::shared_ptr<hid_t> attr_space(new hid_t(id), hemelb::io::hdf5::H5DataSpaceDeleter);
              HEMELB_HDF5_CALL(H5Sget_simple_extent_ndims, rank, *attr_space);
              CPPUNIT_ASSERT_MESSAGE("step offset dataspace is not scalar",
                                     rank == 0);

              // The time dataset doesn't exist within the group
              id = H5Dopen(*group, "time", H5P_DEFAULT);
              CPPUNIT_ASSERT_MESSAGE("time dataset exists", id < 0);
              if (id >= 0)
                H5Dclose(id);
            }

            void testAddSample() {

            }

          private:
            hid_t file, location;
        };

        CPPUNIT_TEST_SUITE_REGISTRATION(TimeDataTests);

      }
    }
  }
}

#endif  // HEMELB_UNITTESTS_IO_H5MD_TIMEDATATESTS_H
