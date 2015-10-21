//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_IO_HDF5_HTRI_H
#define HEMELB_IO_HDF5_HTRI_H

#include <hdf5.h>

namespace hemelb
{
  namespace io
  {
    namespace hdf5
    {

      class HTRI
      {
        public:
          HTRI(const htri_t & tri = 0) :
              tri(tri)
          {
            if (tri < 0)
            {
              throw hdf5::H5Error();
            }
          }

          HTRI & operator=(htri_t & t)
          {
            if (tri != t)
            {
              if (t < 0)
              {
                throw hdf5::H5Error();
              }
              tri = t;
            }
            return *this;
          }

          operator bool() const
          {
            return tri;
          }

        private:
          htri_t tri;
      };

    }
  }
}

#endif  /* HEMELB_IO_HDF5_HTRI_H */
