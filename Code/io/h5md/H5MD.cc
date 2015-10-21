//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#include "io/h5md/H5MD.h"
#include "io/h5md/H5MDError.h"
#include "io/hdf5/H5Types.h"
#include "io/hdf5/H5Error.h"
#include "io/hdf5/HID.h"
#include "io/hdf5/HTRI.h"

#include <hdf5.h>

#define QUOTEx(x) #x
#define QUOTE(x) QUOTEx(x)

namespace hemelb
{
  namespace io
  {
    namespace h5md
    {
      const std::string H5MD::CREATOR_NAME = "hemelb";

      H5MD H5MD::Create(hid_t location, const std::string & creator_name,
                        const std::string & creator_version, const std::string & author_name,
                        const std::string & author_email)
      {
        // Check that a valid identifier has been passed in
        hdf5::HTRI valid = H5Iis_valid(location);
        if (!valid)
        {
          throw H5MDError("Invalid HDF5 identifier passed for location");
        }

        // Check that the identifier corresponds to an HDF5 file or group
        H5I_type_t id_type;
        if ((id_type = H5Iget_type(location)) < 0)
        {
          throw hdf5::H5Error();
        }
        if (id_type != H5I_FILE && id_type != H5I_GROUP)
        {
          throw H5MDError("H5MD structures can only be created within an HDF5 "
                          "file or group");
        }

        // Create an H5MD metadata group
        hdf5::HID<> h5md(H5Gcreate(location, "h5md", H5P_DEFAULT, H5P_DEFAULT,
        H5P_DEFAULT));

        // Attach a version attribute to the H5MD group
        const hsize_t dims[] = { 2 };
        h5md.AttachAttribute("version", VERSION, 1, dims);

        // Create a creator subgroup
        hdf5::HID<> creator(H5Gcreate(h5md, "creator", H5P_DEFAULT, H5P_DEFAULT,
        H5P_DEFAULT));

        // Attach creator name and version attributes to the creator group
        creator.AttachAttribute("name", CREATOR_NAME);
        creator.AttachAttribute("version", QUOTE(HEMELB_REVISION_NUMBER));

        // Create an author subgroup
        hdf5::HID<> author(H5Gcreate(h5md, "author", H5P_DEFAULT, H5P_DEFAULT,
        H5P_DEFAULT));

        // Attach author name and email attributes to the author group
        author.AttachAttribute("name", author_name);
        if (author_email != "")
        {
          author.AttachAttribute("email", author_email);
        }

        return H5MD(location);
      }

      H5MD H5MD::Open(hid_t location)
      {
        if (!IsH5MD(location))
        {
          throw H5MDError("HDF5 location does not contain H5MD structures");
        }
        return H5MD(location);
      }

      bool H5MD::IsH5MD(hid_t location)
      {
        // Check that a valid identifier has been passed in
        hdf5::HTRI valid = H5Iis_valid(location);
        if (!valid)
        {
          throw H5MDError("Invalid HDF5 identifier passed for location");
        }

        // Check that the identifier corresponds to an HDF5 file or group
        H5I_type_t id_type;
        if ((id_type = H5Iget_type(location)) < 0)
        {
          throw hdf5::H5Error();
        }
        if (id_type != H5I_FILE && id_type != H5I_GROUP)
        {
          return false;
        }

        // Check that an object exists with the name "h5md" and it is a group
        hdf5::HTRI exists = H5Oexists_by_name(location, "h5md", H5P_DEFAULT);
        if (!exists)
          return false;

        H5O_info_t info;
        if (H5Oget_info_by_name(location, "h5md", &info, H5P_DEFAULT) < 0)
        {
          throw hdf5::H5Error();
        }
        if (info.type != H5O_TYPE_GROUP)
        {
          return false;
        }

        // Open the H5MD group
        hdf5::HID<> h5md(H5Gopen(location, "h5md", H5P_DEFAULT));

        // Check that the H5MD object has a "version" attribute with a rank 1,
        // length 2 integral dataspace and that the file version is less than or
        // equal to the programme version
        const hsize_t dims[] = { 2 };
        int version[2];
        if (!h5md.CheckAttribute("version", version, 1, dims) ||
            version[0] != VERSION[0] || version[1] > VERSION[1])
        {
          return false;
        }

        // Check that there is a creator group within the h5md group with name
        // and version attributes and that the name is "hemelb"
        exists = H5Oexists_by_name(h5md, "creator", H5P_DEFAULT);
        if (!exists)
        {
          return false;
        }

        if (H5Oget_info_by_name(h5md, "creator", &info, H5P_DEFAULT) < 0)
        {
          throw hdf5::H5Error();
        }
        if (info.type != H5O_TYPE_GROUP)
        {
          return false;
        }

        hdf5::HID<> creator(H5Gopen(h5md, "creator", H5P_DEFAULT));
        std::string name;
        if (!creator.CheckAttribute("name", name) || name != "hemelb" ||
            !creator.CheckAttribute("version"))
        {
          return false;
        }

        // Check that there is an author group within the h5md group with a name
        // attribute
        exists = H5Oexists_by_name(h5md, "author", H5P_DEFAULT);
        if (!exists)
        {
          return false;
        }

        if (H5Oget_info_by_name(h5md, "author", &info, H5P_DEFAULT) < 0)
        {
          throw hdf5::H5Error();
        }
        if (info.type != H5O_TYPE_GROUP)
        {
          return false;
        }

        hdf5::HID<> author(H5Gopen(h5md, "author", H5P_DEFAULT));
        if (!creator.CheckAttribute("name"))
        {
          return false;
        }

        return true;
      }

    }
  }
}
