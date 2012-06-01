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
         * Reads the next bit of data from the file. Returns true if the site was read
         * successfully.
         * @return
         */
        bool ReadNext();

        /**
         * Returns the coordinates of the site.
         * @return
         */
        util::Vector3D<site_t> GetPosition() const;

        /**
         * Returns the pressure at the site.
         * @return
         */
        float GetPressure() const;

        /**
         * Returns the velocity at the site.
         * @return
         */
        util::Vector3D<float> GetVelocity() const;

        /**
         * Returns the shear stress at the site.
         * @return
         */
        float GetShearStress() const;

        /**
         * Returns the Von Mises stress at the site.
         * @return
         */
        float GetVonMisesStress() const;

        /**
         * Returns true if there's a boundary site at given location
         *
         * @param location coordinates of interest
         * @return whether there is a boundary site at location
         */
        bool IsBoundarySite(const util::Vector3D<site_t>& location) const;


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

        /**
         * Variables read in.
         */
        util::Vector3D<site_t> readPosition;
        float readPressure, readShearStress, readVonMisesStress;
        util::Vector3D<float> readVelocity;
    };
  }
}

#endif /* HEMELB_EXTRACTION_SNAPSHOTPARSER_H*/
