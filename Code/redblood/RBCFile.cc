//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#include "redblood/RBCFile.h"
#include "reporting/BuildInfo.h"

namespace hemelb
{
  namespace redblood
  {

    void RBCFile::createH5MD(hid_t file, const std::string & name)
    {
      herr_t error;
      hid_t space, attr;
      hsize_t dims;

      // Create the h5md root group (the only mandatory group in the H5MD spec)
      hid_t h5md;
      if ((h5md = H5Gcreate(file, "h5md", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)) < 0)
        throw;

      // Create a dataspace for the H5MD version attribute
      dims = 2;
      if ((space = H5Screate_simple(1, &dims, nullptr)) < 0)
        throw;

      // Create the H5MD version attribute on the h5md group
      if ((attr = H5Acreate(h5md, "version", H5T_INTEGER, space, H5P_DEFAULT, H5P_DEFAULT)) < 0)
        throw;

      // Write the H5MD version attribute
      if ((error = H5Awrite(attr, H5T_INTEGER, H5MD_VERSION)) < 0)
        throw;

      // Free the H5MD version attribute
      if ((error = H5Aclose(attr)) < 0)
        throw;

      // Free the H5MD version attribute dataspace
      if ((error = H5Sclose(space)) < 0)
        throw;

      // Create the author group as a child of the h5md root group
      hid_t author;
      if ((author = H5Gcreate(h5md, "author", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)) < 0)
        throw;

      // Create a dataspace for the author name attribute
      dims = name.length();
      if ((space = H5Screate_simple(1, &dims, nullptr)) < 0)
        throw;

      // Create the author name attribute on the author group
      if ((attr = H5Acreate(author, "name", H5T_C_S1, space, H5P_DEFAULT, H5P_DEFAULT)) < 0)
        throw;

      // Write the author name attribute
      if ((error = H5Awrite(attr, H5T_C_S1, name.c_str())) < 0)
        throw;

      // Free the author name attribute
      if ((error = H5Aclose(attr)) < 0)
        throw;

      // Free the author name attribute dataspace
      if ((error = H5Sclose(space)) < 0)
        throw;

      // Close the author group
      if ((error = H5Gclose(author)) < 0)
        throw;

      // Create the creator group as a child of the h5md root group
      hid_t creator;
      if ((creator = H5Gcreate(h5md, "creator", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)) < 0)
        throw;

      // Create a dataspace for the creator name attribute
      dims = CREATOR.length();
      if ((space = H5Screate_simple(1, &dims, nullptr)) < 0)
        throw;

      // Create the creator name attribute on the author group
      if ((attr = H5Acreate(creator, "name", H5T_C_S1, space, H5P_DEFAULT, H5P_DEFAULT)) < 0)
        throw;

      // Write the creator name attribute
      if ((error = H5Awrite(attr, H5T_C_S1, CREATOR.c_str())) < 0)
        throw;

      // Free the creator name attribute
      if ((error = H5Aclose(attr)) < 0)
        throw;

      // Free the creator name attribute dataspace
      if ((error = H5Sclose(space)) < 0)
        throw;

      // Create a dataspace for the creator version attribute
      dims = hemelb::reporting::mercurial_revision_number.length();
      if ((space = H5Screate_simple(1, &dims, nullptr)) < 0)
        throw;

      // Create the creator version attribute on the author group
      if ((attr = H5Acreate(creator, "version", H5T_C_S1, space, H5P_DEFAULT, H5P_DEFAULT)) < 0)
        throw;

      // Write the creator version attribute
      if ((error = H5Awrite(attr, H5T_C_S1, hemelb::reporting::mercurial_revision_number.c_str())) < 0)
        throw;

      // Free the creator version attribute
      if ((error = H5Aclose(attr)) < 0)
        throw;

      // Free the creator version attribute dataspace
      if ((error = H5Sclose(space)) < 0)
        throw;

      // Close the author group
      if ((error = H5Gclose(creator)) < 0)
        throw;

      // Close the H5MD root group
      if ((error = H5Gclose(h5md)) < 0)
        throw;
    }

    void RBCFile::checkH5MD(hid_t file)
    {

    }

    std::string RBCFile::readAuthor(hid_t file)
    {

    }

    void RBCFile::write(LatticeTimeStep timestep, PhysicalTime time,
                        const CellContainer & cells)
    {

    }

  }
}
