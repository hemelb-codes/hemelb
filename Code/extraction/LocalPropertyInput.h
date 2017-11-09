
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_EXTRACTION_LOCALPROPERTYINPUT_H
#define HEMELB_EXTRACTION_LOCALPROPERTYINPUT_H

#include "extraction/IterableDataSource.h"
#include "extraction/InputField.h"
//#include "extraction/PropertyOutputFile.h"
#include "net/mpi.h"
#include "net/MpiFile.h"
#include "net/IOCommunicator.h"

namespace hemelb
{
  namespace net
  {
    class IOCommunicator;
  }
  namespace extraction
  {
    /**
     * Stores sufficient information to output property information from this core.
     */
    class LocalPropertyInput
    {
      public:
        /**
         * Initialises a LocalPropertyInput. Required so we can use const reference types.
         */
        LocalPropertyInput(const std::string dataFilePath, const net::IOCommunicator& ioComms);

        /**
         * Tidies up the LocalPropertyInput (close files etc).
         * @return
         */
        ~LocalPropertyInput();

        void LoadDistribution();

	void CheckPreamble();

	void ReadHeaderInfo();

      private:
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

	uint64_t numberOfSites;

        std::vector<InputField> fields;
    };
  }
}

#endif /* HEMELB_EXTRACTION_LOCALPROPERTYOUTPUT_H */
