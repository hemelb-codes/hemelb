#ifndef HEMELB_EXTRACTION_LOCALPROPERTYOUTPUT_H
#define HEMELB_EXTRACTION_LOCALPROPERTYOUTPUT_H

#include "extraction/IterableDataSource.h"
#include "extraction/PropertyOutputFile.h"
#include "mpiInclude.h"

namespace hemelb
{
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
        LocalPropertyOutput(IterableDataSource& dataSource, const PropertyOutputFile& outputSpec);

        /**
         * Tidies up the LocalPropertyOutput (close files etc).
         * @return
         */
        ~LocalPropertyOutput();

        /**
         * Write this core's section of the data file. Only writes if appropriate for the current
         * iteration number
         */
        void Write(unsigned long iterationNumber);

      private:
        /**
         * The MPI file to write into.
         */
        MPI_File outputFile;

        /**
         * The data source to use for file output.
         */
        IterableDataSource& dataSource;

        /**
         * PropertyOutputFile spec.
         */
        const PropertyOutputFile& outputSpec;

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
        char* buffer;
    };
  }
}

#endif /* HEMELB_EXTRACTION_LOCALPROPERTYOUTPUT_H */
