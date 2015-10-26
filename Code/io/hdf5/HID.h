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

      // Destructor functor for HDF5 file handles
      struct H5Fclose
      {
          herr_t operator()(hid_t id)
          {
            return ::H5Fclose(id);
          }
      };

      // Destructor functor for HDF5 group handles
      struct H5Gclose
      {
          herr_t operator()(hid_t id)
          {
            return ::H5Gclose(id);
          }
      };

      // Destructor functor for HDF5 data space handles
      struct H5Sclose
      {
          herr_t operator()(hid_t id)
          {
            return ::H5Sclose(id);
          }
      };

      // Destructor functor for HDF5 dataset handles
      struct H5Dclose
      {
          herr_t operator()(hid_t id)
          {
            return ::H5Dclose(id);
          }
      };

      // Destructor functor for HDF5 attribute handles
      struct H5Aclose
      {
          herr_t operator()(hid_t id)
          {
            return ::H5Aclose(id);
          }
      };

      // Destructor functor for HDF5 property list handles
      struct H5Pclose
      {
          hid_t operator()(hid_t id)
          {
            return ::H5Pclose(id);
          }
      };

      // Destructor functor for HDF5 datatype handles
      struct H5Tclose
      {
          hid_t operator()(hid_t id)
          {
            return ::H5Tclose(id);
          }
      };

      /**
       * Similar to C++11's std::unique_ptr but for HDF5's hid_t type.
       */
      template<class Destructor>
      class HID
      {
        public:

          explicit HID(const hid_t & id, Destructor destructor = Destructor()) :
              id(id), destructor(destructor)
          {
            if (id < 0)
            {
              throw H5Error();
            }
          }

          HID(const HID &) = delete;

          HID(HID && other) :
              id(other.id), destructor(destructor)
          {
            other.id = 0;
          }

          HID & operator=(const HID &) = delete;

          HID & operator=(HID && other)
          {
            if (id != other.id || destructor != other.destructor)
            {
              id = other.id;
              destructor = other.destructor;
              other.id = 0;
            }
            return *this;
          }

          ~HID()
          {
            if (id != 0)
            {
              destructor(id);
            }
          }

          herr_t destroy()
          {
            herr_t error = destructor(id);
            if (error == 0)
            {
              id = 0;
            }
            return error;
          }

          operator hid_t() const
          {
            return id;
          }

        private:
          hid_t id;
          Destructor destructor;
      };

      template<class T>
      void AttachAttribute(hid_t location, const std::string & name, const T * value, int rank,
                           const hsize_t * dims, const hsize_t * maxdims = nullptr, hid_t acpl = H5P_DEFAULT, hid_t aapl = H5P_DEFAULT)
      {
        HID<> space(::H5Screate_simple(rank, dims, maxdims));
        HID<> attribute(::H5Acreate(location, name.c_str(), Types<T>::NATIVE, space, acpl, aapl));
        if (::H5Awrite(attribute, Types<T>::NATIVE, value) < 0)
        {
          throw H5Error();
        }
      }

      void AttachAttribute(hid_t location, const std::string & name, const std::string & value,
                           hid_t acpl = H5P_DEFAULT, hid_t aapl = H5P_DEFAULT)
      {
        HID<> space(::H5Screate(::H5S_SCALAR));
        HID<> type(::H5Tcopy(Types<std::string>::NATIVE));
        if (H5Tset_size(type, std::strlen(value.c_str()) + 1) < 0)
        {
          throw H5Error();
        }
        HID<> attribute(::H5Acreate(location, name.c_str(), type, space, acpl, aapl));
        if (::H5Awrite(attribute, type, value.c_str()) < 0)
        {
          throw H5Error();
        }
      }

      template<class T>
      bool CheckAttribute(hid_t location, const std::string & name, int rank, const hsize_t * dims,
                          hid_t acpl = H5P_DEFAULT, hid_t aapl = H5P_DEFAULT)
      {
        return CheckAttribute<T>(location, name, nullptr, rank, dims, acpl, aapl);
      }

      template<class T>
      bool CheckAttribute(hid_t location, const std::string & name, T * value, int rank,
                          const hsize_t * dims, hid_t acpl = H5P_DEFAULT, hid_t aapl = H5P_DEFAULT)
      {
        // Check the attribute exists
        HTRI exists(::H5Aexists(location, name.c_str()));
        if (!exists)
        {
          return false;
        }

        // Open the attribute and its dataspace
        HID<> attr(::H5Aopen(location, name.c_str(), aapl));
        HID<> space(::H5Aget_space(attr));

        // Check the ranks match
        int actual_rank = ::H5Sget_simple_extent_ndims(space);
        if (actual_rank < 0)
        {
          throw H5Error();
        }
        if (actual_rank != rank)
        {
          return false;
        }

        // Check the dimensions match
        hsize_t actual_dims[rank];
        if (::H5Sget_simple_extent_dims(space, actual_dims, nullptr) < 0)
        {
          throw H5Error();
        }
        if (!std::equal(dims, dims + rank, actual_dims))
        {
          return false;
        }

        // Check the type classes match (the types themselves may not be
        // equal on different platforms but HDF5 is designed to manage this)
        HID<> type(::H5Aget_type(attr));
        if (::H5Tget_class(type) != ::H5Tget_class(Types<T>::NATIVE))
        {
          return false;
        }

        // If the value of the attribute was requested, read it
        if (value != nullptr)
        {
          herr_t error;
          if ((error = ::H5Aread(attr, Types<T>::NATIVE, value)) < 0)
          {
            throw H5Error();
          }
        }

        return true;
      }

      bool CheckAttribute(hid_t location, const std::string & name, std::string & value,
                          hid_t acpl = H5P_DEFAULT, hid_t aapl = H5P_DEFAULT)
      {
        // Check the attribute exists
        HTRI exists(::H5Aexists(location, name.c_str()));
        if (!exists)
        {
          return false;
        }

        // Open the attribute and its dataspace
        HID<> attr(::H5Aopen(location, name.c_str(), aapl));
        HID<> space(::H5Aget_space(attr));

        // Check the ranks match (rank is 0 for strings)
        int rank = ::H5Sget_simple_extent_ndims(space);
        if (rank < 0)
        {
          throw H5Error();
        }
        if (rank != 0)
        {
          return false;
        }

        // Check the type classes match (the types themselves may not be
        // equal on different platforms but HDF5 is designed to manage this)
        HID<> type(::H5Aget_type(attr));
        if (::H5Tget_class(type) != ::H5Tget_class(Types<std::string>::NATIVE))
        {
          return false;
        }

        // Get the size of the string
        std::size_t size;
        if ( (size = ::H5Tget_size(type)) == 0)
        {
          throw H5Error();
        }

        // Read the value of the attribute
        char array[size];
        herr_t error;
        if ( (error = ::H5Aread(attr, Types<std::string>::NATIVE, array)) < 0)
        {
          throw H5Error();
        }
        value.assign(array);

        return true;
      }

      bool CheckAttribute(hid_t location, const std::string & name, hid_t acpl = H5P_DEFAULT,
                          hid_t aapl = H5P_DEFAULT)
      {
        // Check the attribute exists
        HTRI exists = ::H5Aexists(location, name.c_str());
        if (!exists)
        {
          return false;
        }

        // Open the attribute
        HID<> attribute(::H5Aopen(location, name.c_str(), aapl));

        // Check the type classes match (the types themselves may not be
        // equal on different platforms but HDF5 is designed to manage this)
        HID<> type(::H5Aget_type(attribute));
        if (::H5Tget_class(type) != ::H5Tget_class(Types<std::string>::NATIVE))
        {
          return false;
        }

        return true;
      }

    }
  }
}

#endif  /* HEMELB_IO_HDF5_HID_H */
