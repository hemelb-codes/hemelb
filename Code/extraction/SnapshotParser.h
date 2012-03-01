#ifndef HEMELB_EXTRACTION_SNAPSHOTPARSER_H
#define HEMELB_EXTRACTION_SNAPSHOTPARSER_H

#include "io/writers/xdr/XdrFileReader.h"
#include "extraction/IterableDataSource.h"

namespace hemelb
{
  namespace extraction
  {
    class SnapshotParser : public IterableDataSource
    {
      public:
        /**
         * Constructor opens the file, begins reading the header.
         * @param filename
         */
        SnapshotParser(std::string filename);

        /**
         * Destructor ensures file is closed.
         */
        ~SnapshotParser();

        /**
         * Reads the next bit of data from the file, returning the position of the site,
         * and the pressure, velocity and stress there. Returns true if the site was read
         * successfully.
         * @param position
         * @param pressure
         * @param velocity
         * @param stress
         * @return
         */
        bool ReadNext(util::Vector3D<float>& position, float& pressure, util::Vector3D<float>& velocity, float& stress);

      private:
        /**
         * File handle to open snapshot.
         */
        FILE* file;
        /**
         * XDR reader for decoding the file into memory.
         */
        io::writers::xdr::XdrFileReader reader;

        /**
         * The length in real-world units of each lattice unit.
         */
        double voxelSize;
        /**
         * The real-world origin of the simulation coordinate system.
         */
        util::Vector3D<double> origin;
        /**
         * The coordinates of the lower-left-nearest site in the simulation.
         */
        util::Vector3D<int> minimalSite;
        /**
         * The number of fluid sites.
         */
        int fluidSiteCount;
        /**
         * Count of the number of sites read so far.
         */
        int sitesRead;
    };
  }
}

#endif /* HEMELB_EXTRACTION_SNAPSHOTPARSER_H*/
