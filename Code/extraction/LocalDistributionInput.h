
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_EXTRACTION_LOCALDISTRIBUTIONINPUT_H
#define HEMELB_EXTRACTION_LOCALDISTRIBUTIONINPUT_H

#include "extraction/IterableDataSource.h"
#include "extraction/InputField.h"
#include "io/writers/xdr/XdrMemReader.h"
#include "lb/lattices/Lattices.h"
#include "net/mpi.h"
#include "net/MpiFile.h"
#include "net/IOCommunicator.h"

namespace hemelb
{
  namespace net
  {
    class IOCommunicator;
  }
  namespace geometry
  {
    class LatticeData;
  }
  namespace extraction
  {
    /**
     * Stores sufficient information to output property information from this core.
     */
    class LocalDistributionInput
    {
      public:
        /**
         * Initialises a LocalDistributionInput. Required so we can use const reference types.
         */
        LocalDistributionInput(const std::string dataFilePath, const net::IOCommunicator& ioComms);

        /**
         * Tidies up the LocalDistributionInput (close files etc).
         * @return
         */
        ~LocalDistributionInput();

        void LoadDistribution(geometry::LatticeData* latDat);

      private:
	typedef hemelb::lb::lattices:: HEMELB_LATTICE LatticeType;

	void CheckPreamble();

	void ReadHeaderInfo();

        const net::IOCommunicator& comms;

        // The rank which reads in the header information.
        static const proc_t HEADER_READING_RANK = 0;

        /**
         * The path to the file to read from.
         */
        const std::string filePath;

        /**
         * Where to begin reading from the file.
         */
        //uint64_t localDataOffsetIntoFile;

        /**
         * The length, in bytes, of the local read.
         */
        //uint64_t readLength;

        /**
         * The length, in bytes, of the total read length;
         */
        //uint64_t allCoresReadLength;

        /**
         * Buffer to write into before writing to disk.
         */
        //std::vector<char> buffer;

        /**
         * Type of written values
         */
        //typedef float WrittenDataType;

        // File accessed to read in the distributions.
        net::MpiFile inputFile;

        // File accessed to read in the offsets of the distributions.
        net::MpiFile offsetFile;

	uint64_t numberOfSites;

        std::vector<InputField> fields;

        uint64_t thisOffset;

        uint64_t nextOffset;

	uint32_t lengthOfSegment;

	std::vector<char>* dataBufferPtr;

	io::writers::xdr::XdrReader* dataReaderPtr;
    };
  }
}

#endif /* HEMELB_EXTRACTION_LOCALDISTRIBUTIONINPUT_H */
