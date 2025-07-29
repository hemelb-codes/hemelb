// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

// abstract class LocalPropertyOutput
// derived classes are XdrPropertyOutput and Hdf5PropertyOutput

#ifndef HEMELB_EXTRACTION_LOCALPROPERTYOUTPUT_H
#define HEMELB_EXTRACTION_LOCALPROPERTYOUTPUT_H

#include "extraction/IterableDataSource.h"
#include "extraction/PropertyOutputFile.h"
#include "lb/Lattices.h"
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
    // Stores sufficient information to output property information
    // from this core.
    class LocalPropertyOutput
    {
    public:
      // Initialises a LocalPropertyOutput. Required so we can use
      // const reference types. Collective on the communicator.
      LocalPropertyOutput(IterableDataSource& dataSource, const PropertyOutputFile& outputSpec,
			  const net::IOCommunicator& ioComms);
      // a virtual destructor function for abstract class
      virtual ~LocalPropertyOutput() = default;

      // True if this property output should be written on the current iteration.
      bool ShouldWrite(unsigned long timestepNumber) const;

      // Returns the property output file object to be written.
      const PropertyOutputFile& GetOutputSpec() const;

      // Write this core's section of the data file. Only writes if
      // appropriate for the current iteration number
      // only work as interface specification
      virtual void Write(unsigned long timestepNumber, unsigned long totalSteps)=0;

    protected:
      // How many sites does this MPI process write?
      std::uint64_t CountWrittenSitesOnRank();
      // Returns the number of items written for the field.
      unsigned GetFieldLength(source::Type) const;
      // Our communicator
      const net::IOCommunicator& comms;
      // The data source to use for file output.
      IterableDataSource& dataSource;
      // PropertyOutputFile spec.
      PropertyOutputFile outputSpec;

      std::uint64_t local_site_count;
      std::uint64_t global_site_count;

    private:
    };
  }
}

#endif // HEMELB_EXTRACTION_LOCALPROPERTYOUTPUT_H
