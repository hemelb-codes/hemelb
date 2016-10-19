
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_COMM_MPIFILE_H
#define HEMELB_COMM_MPIFILE_H

#include <memory>
#include "comm/MpiError.h"
#include "comm/MpiCommunicator.h"

namespace hemelb
{
  namespace comm
  {
    class MpiFile
    {
      public:
        MpiFile();

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

        void SetView(MPI_Offset disp, MPI_Datatype etype, MPI_Datatype filetype,
		     const std::string& datarep, const MPI_Info info = MPI_INFO_NULL);

	MPI_Offset GetSize() const;

        const MpiCommunicator& GetCommunicator() const;

        template<typename T>
        void Read(std::vector<T>& buffer, MPI_Status* stat = MPI_STATUS_IGNORE);
        template<typename T>
        void ReadAt(MPI_Offset offset, std::vector<T>& buffer, MPI_Status* stat = MPI_STATUS_IGNORE);

        template<typename T>
        void Write(const std::vector<T>& buffer, MPI_Status* stat = MPI_STATUS_IGNORE);
        template<typename T>
        void WriteAt(MPI_Offset offset, const std::vector<T>& buffer, MPI_Status* stat = MPI_STATUS_IGNORE);
      private:
	friend class MpiCommunicator;
	
        MpiFile(const MpiCommunicator* parentComm, MPI_File fh);

        const MpiCommunicator* comm;
        std::shared_ptr<MPI_File> filePtr;
    };

  }
}

#include "comm/MpiFile.hpp"

#endif /* HEMELB_COMM_MPIFILE_H */
