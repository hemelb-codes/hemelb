// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

// LocalPropertyHdf5Output.h
// In extraction/LocalPropertyHdf5Output.h
// one hdf5 file, group, separate datasets, and XDMF

#pragma once

#ifdef USE_HDF5

#include "extraction/LocalPropertyOutput.h"
#include <hdf5.h>
#include <string>
#include <vector>

namespace hemelb
{
  namespace extraction
  {
    class LocalPropertyHdf5Output final : public LocalPropertyOutput
    {
    public:
      LocalPropertyHdf5Output(IterableDataSource& dataSource,
                              const PropertyOutputFile& outputSpec,
                              const net::IOCommunicator& ioComms);
      
      ~LocalPropertyHdf5Output();
      
      void Write(unsigned long timestepNumber, unsigned long totalSteps) override;

    private:
      // Writes the XDMF file for the HDF5 output.
      void WriteXDMFFile();
      
      MPI_Comm mpi_comm;
    
      hid_t file_id = -1;
      
      // vector for written timesteps
      std::vector<unsigned long> written_timesteps;
    };
  }
}
#endif // USE_HDF5