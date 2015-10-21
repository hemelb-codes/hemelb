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

#include <exception>
#include <string>

namespace hemelb
{
  namespace io
  {

    namespace hdf5
    {

      /**
       * Indicate an error to do with HDF5.
       */
      class H5Error : public std::exception
      {
        public:
          H5Error();

          virtual const char * what() const noexcept
          {
            return wat.c_str();
          }

        private:
          std::string wat;
      };

    }
  }
}

#endif // HEMELB_IO_HDF5_H5ERROR_H
