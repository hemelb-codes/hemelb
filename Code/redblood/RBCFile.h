//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_REDBLOOD_RBCFILE_H
#define HEMELB_REDBLOOD_RBCFILE_H

#include <hdf5.h>
#include <mpi.h>

#include "redblood/Cell.h"

namespace hemelb
{
  namespace redblood
  {

    class RBCFile
    {
      public:
        static const std::string CREATOR = "HemeLB";
        static const int H5MD_VERSION[2] = { 1, 0 };

        /**
         * Creates a new RBC file and opens it for read/write access.
         *
         * @param filename the name of the file to create.
         * @param trunc  if true and the file already exists, truncate the
         *   existing contents and open with read/write access.  If false fail
         *   if the file already exists.
         * @throw if the file already exists.
         * @return a shared pointer to an RBCFile object.
         */
        static std::shared_ptr<RBCFile> create(const std::string & filename,
                                               const std::string & author,
                                               bool trunc = false)
        {
          // Create a new HDF5 file
          hid_t file_id = H5Fcreate(filename.c_str(),
                                    (trunc) ? H5F_ACC_TRUNC : H5F_ACC_EXCL,
                                    H5P_DEFAULT, H5P_DEFAULT);

          // If creation failed throw an exception
          if (file_id < 0)
            throw;

          // Create H5MD standard groups/datasets
          createH5MD(file_id, author);

          // Return a shared pointer to a newly created RBCFile object
          return std::make_shared<RBCFile>(file_id);
        }

        /**
         * Creates a new RBC file and opens it for read/write access using MPI
         * IO.  This is a collective operation.
         *
         * @param filename the name of the file to create.
         * @param comm the MPI IO communicator.
         * @param info the MPI Info object to be used to open the file.
         * @param trunc  if true and the file already exists, truncate the
         *   existing contents and open with read/write access.  If false fail
         *   if the file already exists.
         * @throw if the file already exists.
         * @return a shared pointer to an RBCFile object.
         */
        static std::shared_ptr<RBCFile> create(const std::string & filename,
                                               const std::string & author,
                                               const MPI_Comm & comm,
                                               const MPI_Info & info,
                                               bool trunc = false)
        {
          // TODO: This method signature will probably have to be changed when
          // MPI is implemented for RBCs to take a hemelb::net::MpiCommunicator
          // or whatever the chosen implementation is
          herr_t status;
          hid_t access_plist = H5Pcreate(H5P_FILE_ACCESS);
          if (access_plist < 0)
            throw;

          if ((status = H5Pset_fapl_mpio(access_plist, comm, info)) < 0)
            throw;

          // Create a new HDF5 file
          hid_t file_id = H5Fcreate(filename.c_str(),
                                    (trunc) ? H5F_ACC_TRUNC : H5F_ACC_EXCL,
                                    H5P_DEFAULT, access_plist);

          // If creation failed throw an exception
          if (file_id < 0)
            throw;

          // Create H5MD standard groups/datasets
          createH5MD(file_id, author);

          // Return a shared pointer to a newly created RBCFile object
          return std::make_shared<RBCFile>(file_id);
        }

        /**
         * Opens an existing RBC file.
         *
         * @param filename the name of the file to open.
         * @param readonly whether to open the file with readonly access (default)
         *   or read/write.
         * @throw if the file does not exist.
         * @return a shared pointer to an RBCFile object.
         */
        static std::shared_ptr<RBCFile> open(const std::string & filename,
                                             bool readonly = true)
        {
          // Open the HDF5 file
          hid_t file_id = H5Fopen(filename.c_str(),
                                  (readonly) ? H5F_ACC_RDONLY : H5F_ACC_RDWR,
                                  H5P_DEFAULT);

          // If opening failed throw an exception
          if (file_id < 0)
            throw;

          // Check for the H5MD standard groups/datasets
          checkH5MD(file_id);

          // Return a shared pointer to a newly created RBCFile object
          return std::make_shared<RBCFile>(file_id);
        }

        /**
         * Opens an existing RBC file using MPI IO.  This is a collective
         * operation.
         *
         * @param filename the name of the file to open.
         * @param comm the MPI IO communicator.
         * @param info the MPI Info object to be used to open the file.
         * @param readonly whether to open the file with readonly access (default)
         *   or read/write.
         * @throw if the file does not exist.
         * @return a shared pointer to an RBCFile object.
         */
        static std::shared_ptr<RBCFile> open(const std::string & filename,
                                             const MPI_Comm & comm = MPI_COMM_WORLD,
                                             const MPI_Info & info = MPI_INFO_NULL,
                                             bool readonly = true)
        {
          // TODO: This method signature will probably have to be changed when
          // MPI is implemented for RBCs to take a hemelb::net::MpiCommunicator
          // or whatever the chosen implementation is
          herr_t status;
          hid_t access_plist = H5Pcreate(H5P_FILE_ACCESS);
          if (access_plist < 0)
            throw;

          if ( (status = H5Pset_fapl_mpio(access_plist, comm, info)) < 0)
            throw;

          // Open the HDF5 file
          hid_t file_id = H5Fopen(filename.c_str(),
                                  (readonly) ? H5F_ACC_RDONLY : H5F_ACC_RDWR,
                                  access_plist);

          // If opening failed throw an exception
          if (file_id < 0)
            throw;

          // Check for the H5MD standard groups/datasets
          checkH5MD(file_id);

          // Return a shared pointer to a newly created RBCFile object
          return std::make_shared<RBCFile>(file_id);
        }

        ~RBCFile()
        {
          // Try to close the file and prevent any exception from escaping
          try
          {
            close();
          }
          catch (std::exception & e)
          {
            // TODO: Log the exception
          }
        }

        void write(LatticeTimeStep, PhysicalTime, const CellContainer &);

        void flush()
        {
          herr_t status = H5Fflush(file, H5F_SCOPE_LOCAL);
          if (status < 0)
            throw;
        }

        /**
         * Closes the RBC file.  This is a collective operation if the file
         * was opened using MPI IO.
         */
        void close()
        {
          herr_t status = H5Fclose(file);
          if (status < 0)
            throw;
        }

        void SetWriteBaryCentre(bool b)
        {
          writeBaryCentre = b;
        }
        bool GetWriteBaryCentre()
        {
          return writeBaryCentre;
        }
        void SetWriteMesh(bool b)
        {
          writeMesh = b;
        }
        bool GetWriteMesh()
        {
          return writeMesh;
        }
        void SetWriteScale(bool b)
        {
          writeScale = b;
        }
        bool GetWriteScale()
        {
          return writeScale;
        }

        std::string GetAuthor()
        {
          return author;
        }

      private:

        static void createH5MD(hid_t, const std::string &);
        static void checkH5MD(hid_t);
        std::string readAuthor(hid_t);

        RBCFile(const hid_t & id) :
            file(id), writeBaryCentre(true), writeMesh(false), writeScale(false),
            author(readAuthor(id))
        {
        }

        hid_t file;
        bool writeBaryCentre, writeMesh, writeScale;
        const std::string author;
    };

  }
}

#endif
