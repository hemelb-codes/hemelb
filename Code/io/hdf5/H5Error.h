//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_IO_HDF5_H5ERROR_H
#define HEMELB_IO_HDF5_H5ERROR_H

#include <hdf5.h>
#include "Exception.h"

#include <string>

namespace hemelb
{
  namespace io
  {

    namespace hdf5
    {
      /**
       * Indicate an error to do with HDF5.
       *
       * Will be thrown by the HEMELB_HDF5_CALL macro, as in:
       *   HEMELB_HDF5_CALL(H5Fcreate, file_id, filename, flags, create_plist, access_plist);
       *
       */
      class H5Error : public ::hemelb::Exception
      {
        public:
          H5Error(const std::string &, herr_t, const std::string &, unsigned int);

          const std::string & GetFunction() {
            return function;
          }

          hid_t GetError() {
            return error;
          }

          const std::string & GetFile() {
            return file;
          }

          unsigned int GetLine() {
            return line;
          }

        private:
          const std::string & function;
          const herr_t error;
          const std::string & file;
          const unsigned int line;
      };

      void H5FileDeleter(hid_t *);
      void H5GroupDeleter(hid_t *);
      void H5DataSetDeleter(hid_t *);
      void H5DataSpaceDeleter(hid_t *);
      void H5AttributeDeleter(hid_t *);
      void H5PropertyListDeleter(hid_t *);
    }
  }
}

#define HEMELB_HDF5_CALL(function, res, ...) \
do \
{ \
  if ((res = function(__VA_ARGS__)) < 0) \
    throw ::hemelb::io::hdf5::H5Error(#function, res, __FILE__, __LINE__); \
} while (false)

#endif // HEMELB_IO_HDF5_H5ERROR_H
