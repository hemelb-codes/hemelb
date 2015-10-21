//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_IO_HDF5_HID_H
#define HEMELB_IO_HDF5_HID_H

#include <hdf5.h>

#include "log/Logger.h"

namespace hemelb
{
  namespace io
  {
    namespace hdf5
    {

      struct H5Destroy
      {
          herr_t operator()(hid_t id)
          {
            H5I_type_t type = H5Iget_type(id);
            herr_t error = -1;
            switch (type)
            {
              case H5I_FILE:
                error = H5Fclose(id);
                break;
              case H5I_GROUP:
                error = H5Gclose(id);
                break;
              case H5I_DATATYPE:
                error = H5Tclose(id);
                break;
              case H5I_DATASPACE:
                error = H5Sclose(id);
                break;
              case H5I_DATASET:
                error = H5Dclose(id);
                break;
              case H5I_ATTR:
                error = H5Aclose(id);
                break;
              default:
                hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("Closing HDF5 identifier failed: unknown type");
            }
            return error;
          }
      };

      /**
       * Similar to C++11's std::lock_guard but for HDF5's hid_t type.
       */
      template<class Destructor = H5Destroy>
      class HID
      {
        public:
          explicit HID(const hid_t & id, const Destructor & destroy = Destructor()) :
              id(id), destroy(destroy)
          {
            if (id < 0)
            {
              throw hdf5::H5Error();
            }
          }

          HID(const HID &) = delete;

          ~HID()
          {
            destroy(id);
          }

          operator hid_t() const
          {
            return id;
          }

          template<class T>
          void AttachAttribute(const char * name, const T * value, int rank, const hsize_t * dims,
                               hid_t acpl = H5P_DEFAULT, hid_t aapl = H5P_DEFAULT)
          {
            hdf5::HID<> space(H5Screate_simple(rank, dims, dims));
            hdf5::HID<> attr(H5Acreate(id, name, hdf5::Types<T>::NATIVE, space, acpl, aapl));
            if (H5Awrite(attr, hdf5::Types<T>::NATIVE, value) < 0)
            {
              throw hdf5::H5Error();
            }
          }

          void AttachAttribute(const char * name, const std::string & value, hid_t acpl =
          H5P_DEFAULT,
                               hid_t aapl = H5P_DEFAULT)
          {
            hdf5::HID<> space(H5Screate(H5S_SCALAR));
            hdf5::HID<> type(H5Tcopy(hdf5::Types<std::string>::NATIVE));
            if (H5Tset_size(type, std::strlen(value.c_str()) + 1) < 0)
            {
              throw hdf5::H5Error();
            }
            hdf5::HID<> attr(H5Acreate(id, name, type, space, acpl, aapl));
            if (H5Awrite(attr, type, value.c_str()) < 0)
            {
              throw hdf5::H5Error();
            }
          }

          template<class T>
          bool CheckAttribute(const char * name, T * value, int rank, const hsize_t * dims,
                              hid_t acpl = H5P_DEFAULT, hid_t aapl = H5P_DEFAULT)
          {
            // Check the attribute exists
            hdf5::HID<> exists(H5Aexists(id, name));
            if (!exists)
            {
              return false;
            }

            // Open the attribute and its dataspace
            hdf5::HID<> attr(H5Aopen(id, name, aapl));
            hdf5::HID<> space(H5Aget_space(attr));

            // Check the ranks match
            int actual_rank = H5Sget_simple_extent_ndims(space);
            if (actual_rank < 0)
            {
              throw hdf5::H5Error();
            }
            if (actual_rank != rank)
            {
              return false;
            }

            // Check the dimensions match
            hsize_t actual_dims[rank];
            if (H5Sget_simple_extent_dims(space, actual_dims, nullptr) < 0)
            {
              throw hdf5::H5Error();
            }
            if (!std::equal(dims, dims + rank, actual_dims))
            {
              return false;
            }

            // Check the type classes match (the types themselves may not be
            // equal on different platforms but HDF5 is designed to manage this)
            hdf5::HID<> type(H5Aget_type(attr));
            if (H5Tget_class(type) != H5Tget_class(hdf5::Types<T>::NATIVE))
            {
              return false;
            }

            // If the value of the attribute was requested, read it
            if (value != nullptr)
            {
              herr_t error;
              if ((error = H5Aread(attr, hdf5::Types<T>::NATIVE, value)) < 0)
              {
                throw hdf5::H5Error();
              }
            }

            return true;
          }

          bool CheckAttribute(const char * name, std::string & value, hid_t acpl = H5P_DEFAULT,
                              hid_t aapl = H5P_DEFAULT)
          {
            // Check the attribute exists
            hdf5::HID<> exists(H5Aexists(id, name));
            if (!exists)
            {
              return false;
            }

            // Open the attribute and its dataspace
            hdf5::HID<> attr(H5Aopen(id, name, aapl));
            hdf5::HID<> space(H5Aget_space(attr));

            // Check the ranks match (rank is 0 for strings)
            int rank = H5Sget_simple_extent_ndims(space);
            if (rank < 0)
            {
              throw hdf5::H5Error();
            }
            if (rank != 0)
            {
              return false;
            }

            // Check the type classes match (the types themselves may not be
            // equal on different platforms but HDF5 is designed to manage this)
            hdf5::HID<> type(H5Aget_type(attr));
            if (H5Tget_class(type) != H5Tget_class(hdf5::Types<std::string>::NATIVE))
            {
              return false;
            }

            // Get the size of the string
            std::size_t size;
            if ( (size = H5Tget_size(type)) == 0)
            {
              throw hdf5::H5Error();
            }

            // Read the value of the attribute
            char array[size];
            herr_t error;
            if ( (error = H5Aread(attr, hdf5::Types<std::string>::NATIVE, array)) < 0)
            {
              throw hdf5::H5Error();
            }
            value.assign(array);

            return true;
          }

          bool CheckAttribute(const char * name, hid_t acpl = H5P_DEFAULT, hid_t aapl = H5P_DEFAULT)
          {
            // Check the attribute exists
            hdf5::HID<> exists(H5Aexists(id, name));
            if (!exists)
              return false;

            // Open the attribute
            hdf5::HID<> attr(H5Aopen(id, name, aapl));

            // Check the type classes match (the types themselves may not be
            // equal on different platforms but HDF5 is designed to manage this)
            hdf5::HID<> type(H5Aget_type(attr));
            if (H5Tget_class(type) != H5Tget_class(hdf5::Types<std::string>::NATIVE))
            {
              return false;
            }

            return true;
          }

        private:
          hid_t id;
          Destructor destroy;
      };

    }
  }
}

#endif  /* HEMELB_IO_HDF5_HID_H */
