
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_EXTRACTION_LOCALDISTRIBUTIONINPUT_H
#define HEMELB_EXTRACTION_LOCALDISTRIBUTIONINPUT_H

#include <boost/optional.hpp>

#include "comm/Communicator.h"
#include "extraction/IterableDataSource.h"
#include "extraction/InputField.h"
#include "io/writers/xdr/XdrMemReader.h"
#include "lb/lattices/Lattices.h"

namespace hemelb
{
  namespace comm
  {
    class Communicator;
    class MpiFile;
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
      LocalDistributionInput(const std::string& dataFilePath, comm::Communicator::ConstPtr ioComms);

        /**
         * Tidies up the LocalDistributionInput (close files etc).
         * @return
         */
        ~LocalDistributionInput();

        void LoadDistribution(geometry::LatticeData* latDat, boost::optional<LatticeTimeStep>& initalTime);

      private:
	typedef hemelb::lb::lattices:: HEMELB_LATTICE LatticeType;

        void ReadExtractionHeaders(comm::MpiFile&);
        void ReadOffsets(const std::string&);

        comm::Communicator::ConstPtr comms;

        /**
         * The path to the file to read from.
         */
        std::string filePath;

        InputField distField;
        uint64_t localStart;
        uint64_t localStop;
        uint64_t timestep;
        uint64_t allCoresWriteLength;
    };
  }
}

#endif /* HEMELB_EXTRACTION_LOCALDISTRIBUTIONINPUT_H */
