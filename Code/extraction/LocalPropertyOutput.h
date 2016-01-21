
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_EXTRACTION_LOCALPROPERTYOUTPUT_H
#define HEMELB_EXTRACTION_LOCALPROPERTYOUTPUT_H

#include "extraction/IterableDataSource.h"
#include "extraction/PropertyOutputFile.h"
#include "net/mpi.h"
#include "net/MpiFile.h"

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
    class LocalPropertyOutput
    {
      public:
        /**
         * Initialises a LocalPropertyOutput. Required so we can use const reference types.
         * @param file
         * @param offset
         * @return
         */
        LocalPropertyOutput(IterableDataSource& dataSource, const PropertyOutputFile* outputSpec, const net::IOCommunicator& ioComms);

        /**
         * Tidies up the LocalPropertyOutput (close files etc).
         * @return
         */
        ~LocalPropertyOutput();

        /**
         * True if this property output should be written on the current iteration.
         * @return
         */
        bool ShouldWrite(unsigned long timestepNumber) const;

        /**
         * Returns the property output file object to be written.
         * @return
         */
        const PropertyOutputFile* GetOutputSpec() const;

        /**
         * Write this core's section of the data file. Only writes if appropriate for the current
         * iteration number
         */
        void Write(unsigned long timestepNumber);

      private:
        /**
         * Returns the number of floats written for the field.
         * @param field
         */
        unsigned GetFieldLength(OutputField::FieldType field);

        /**
         * Returns the offset to the field, as it should be written to file.
         * @param field
         * @return
         */
        double GetOffset(OutputField::FieldType field) const;

        const net::IOCommunicator& comms;
        /**
         * The MPI file to write into.
         */
        net::MpiFile outputFile;

        /**
         * The data source to use for file output.
         */
        IterableDataSource& dataSource;

        /**
         * PropertyOutputFile spec.
         */
        const PropertyOutputFile* outputSpec;

        /**
         * Where to begin writing into the file.
         */
        uint64_t localDataOffsetIntoFile;

        /**
         * The length, in bytes, of the local write.
         */
        uint64_t writeLength;

        /**
         * The length, in bytes, of the total write length;
         */
        uint64_t allCoresWriteLength;

        /**
         * Buffer to write into before writing to disk.
         */
        std::vector<char> buffer;

        /**
         * Type of written values
         */
        typedef float WrittenDataType;
    };
  }
}

#endif /* HEMELB_EXTRACTION_LOCALPROPERTYOUTPUT_H */
