//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_NET_MPIFILE_H
#define HEMELB_NET_MPIFILE_H

#include <boost/shared_ptr.hpp>
#include "net/MpiError.h"
#include "net/MpiCommunicator.h"

namespace hemelb
{
  namespace net
  {
    class MpiFile
    {
      public:
        MpiFile();

        /**
         * Opens a file with MPI_File_open. A collective operation on comm.
         * @param comm
         * @param filename
         * @param mode
         * @param info
         * @return
         */
        static MpiFile Open(const MpiCommunicator& comm, const std::string& filename, int mode,
                            const MPI_Info info = MPI_INFO_NULL);

        /**
         * Closes the file with MPI_File_close.
         * Note this is a collective operation.
         */
        void Close();

        /**
         * Allow implicit casts to MPI_file
         * @return The underlying MPI file handle
         */
        operator MPI_File() const;

        void SetView(MPI_Offset disp, MPI_Datatype etype, MPI_Datatype filetype, const std::string& datarep, MPI_Info info);

        const MpiCommunicator& GetCommunicator() const;

        template<typename T>
        void Read(std::vector<T>& buffer, MPI_Status* stat = MPI_STATUS_IGNORE);
        template<typename T>
        void ReadAt(MPI_Offset offset, std::vector<T>& buffer, MPI_Status* stat = MPI_STATUS_IGNORE);

        template<typename T>
        void Write(std::vector<T>& buffer, MPI_Status* stat = MPI_STATUS_IGNORE);
        template<typename T>
        void WriteAt(MPI_Offset offset, std::vector<T>& buffer, MPI_Status* stat = MPI_STATUS_IGNORE);
      protected:
        MpiFile(const MpiCommunicator& parentComm, MPI_File fh);

        const MpiCommunicator* comm;
        boost::shared_ptr<MPI_File> filePtr;
    };

  }
}

#include "net/MpiFile.hpp"

#endif /* HEMELB_NET_MPIFILE_H */
