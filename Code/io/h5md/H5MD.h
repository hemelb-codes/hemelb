//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_IO_H5MD_H5MD_H
#define HEMELB_IO_H5MD_H5MD_H

#include <string>
#include <memory>
#include <hdf5.h>

#include "io/hdf5/H5Error.h"

namespace hemelb
{
  namespace io
  {
    namespace h5md
    {

      class H5MD
      {

        public:

          static constexpr int VERSION[] = { 1, 1 };
          static const std::string CREATOR_NAME;

          /*!
           * Creates an H5MD structure rooted at `location` within an existing
           * HDF5 file.  `location` may be the root of the file or an existing
           * group within the file.
           *
           * \param location         where to create the H5MD structures
           * \param creator_name     the name of the program writing the H5MD
           *                           file
           * \param creator_version  the version string of the program writing
           *                           the H5MD file
           * \param author_name      the name of the person performing the
           *                           simulation
           * \param author_email     the email address of the person performing
           *                           the simulation
           * \throw H5MDError if the location passed in is not a valid HDF5 file
           *                    or group
           * \return an H5MD object.
           */
          static std::shared_ptr<H5MD> Create(hid_t, const std::string &,
                                              const std::string &,
                                              const std::string &,
                                              const std::string & = "");

          /*!
           * Opens an existing H5MD structure rooted at `location`.
           *
           * \param location  the root of the H5MD structure
           * \throw H5MDError if the location passed in does not represent a
           *                    valid H5MD location (see #IsH5MD)
           * \return an H5MD object.
           */
          static std::shared_ptr<H5MD> Open(hid_t);

          /*!
           * Checks whether the location within an HDF5 file contains an H5MD
           * structure.
           *
           * \param location  the root of the H5MD structure
           * \throw H5MDError if the location passed in is not a valid HDF5
           *                    identifier
           * \return true if the location is a valid H5MD structure, false
           *           otherwise
           */
          static bool IsH5MD(hid_t);

        private:
          H5MD(hid_t location) :
              location(location)
          {
          }
          hid_t location;
      };

    }
  }
}

#endif  // HEMELB_IO_H5MD_H5MD_H
