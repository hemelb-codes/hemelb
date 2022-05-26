// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_EXTRACTION_LOCALDISTRIBUTIONINPUT_H
#define HEMELB_EXTRACTION_LOCALDISTRIBUTIONINPUT_H

#include <optional>

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
    // Read this ranks's part of a checkpoint file.
    class LocalDistributionInput
    {
      public:
      // Construct, but don't do any I/O.
      //
      // Take the path string by value since we will move it into a member
      // anyway.
      LocalDistributionInput(std::string dataFilePath, const net::IOCommunicator& ioComms);

      // Open the file and load our part into the LatticeData
      // instance.
      //
      // Time is optional, if not supplied will use the last one in
      // the file.
      //
      // Requires the checkpoint have been saved with exactly the same
      // domain decomposition as currently running.
      void LoadDistribution(geometry::LatticeData* latDat, std::optional<LatticeTimeStep>& initalTime);

    private:

      void ReadExtractionHeaders(net::MpiFile&, const unsigned NUMVECTORS);
      void ReadOffsets(const std::string&);

      const net::IOCommunicator& comms;

      // The path to the file to read from.
      std::string filePath;

      InputField distField;
      uint64_t localStart;
      uint64_t localStop;
      uint64_t timestep;
      uint64_t allCoresWriteLength;
    };
  }
}

#endif // HEMELB_EXTRACTION_LOCALDISTRIBUTIONINPUT_H
