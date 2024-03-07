// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_NET_MPIFILE_H
#define HEMELB_NET_MPIFILE_H

#include <filesystem>
#include <memory>
#include <span>

#include "net/MpiError.h"
#include "net/MpiCommunicator.h"

namespace hemelb::net
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
        static MpiFile Open(const MpiCommunicator& comm, const std::filesystem::path& filename, int mode,
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

        void SetView(MPI_Offset disp, MPI_Datatype etype, MPI_Datatype filetype,
		     const std::string& datarep, const MPI_Info info = MPI_INFO_NULL);

	MPI_Offset GetSize() const;

        const MpiCommunicator& GetCommunicator() const;

        template<typename T, std::size_t N>
        void Read(std::span<T, N> buffer, MPI_Status* stat = MPI_STATUS_IGNORE);
        template<typename T, std::size_t N>
        void ReadAt(MPI_Offset offset, std::span<T, N> buffer,
                    MPI_Status* stat = MPI_STATUS_IGNORE);
        template<typename T, std::size_t N>
        void ReadAtAll(MPI_Offset offset, std::span<T, N> buffer,
                       MPI_Status* stat = MPI_STATUS_IGNORE);

        template<typename T, std::size_t N>
        void Write(std::span<T const, N> buffer, MPI_Status* stat = MPI_STATUS_IGNORE);
        template<typename T, std::size_t N>
        void WriteAt(MPI_Offset offset, std::span<T const, N> buffer, MPI_Status* stat =
                         MPI_STATUS_IGNORE);
    protected:
        MpiFile(const MpiCommunicator& parentComm, MPI_File fh);

        const MpiCommunicator* comm;
        std::shared_ptr<MPI_File> filePtr;
    };

}

#include "net/MpiFile.hpp"

#endif /* HEMELB_NET_MPIFILE_H */
