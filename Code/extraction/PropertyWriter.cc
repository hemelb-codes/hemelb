// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "extraction/PropertyWriter.h"
#include "extraction/XdrPropertyOutput.h"
#include "extraction/LocalPropertyHdf5Output.h"

namespace hemelb
{
  namespace extraction
  {
    PropertyWriter::PropertyWriter(IterableDataSource& dataSource,
                                   const std::vector<PropertyOutputFile>& propertyOutputs,
                                   const net::IOCommunicator& ioComms)
    {
      for (unsigned outputNumber = 0; outputNumber < propertyOutputs.size(); ++outputNumber)
      {
        // Get the current format specification
        const auto& currentOutputSpec = propertyOutputs[outputNumber];
        const std::string& format = currentOutputSpec.format;
        // Create the appropriate output type based on the format.
        if (format == "hdf5")
        {
          #ifdef USE_HDF5
          //create HDF5 output
          localPropertyOutputs.push_back(new LocalPropertyHdf5Output(dataSource,
                                                                     currentOutputSpec,
                                                                     ioComms));
          #else
          throw Exception() << "HDF5 output was requested in the XML configuration, but HDF5 support is not enabled in this build."
                            << "Please recompile with -DUSE_HDF5=ON.";
          #endif
        }
        else if(format == "xdr")
        {
          //create XDR output
          localPropertyOutputs.push_back(new XdrPropertyOutput(dataSource,
                                                               currentOutputSpec,
                                                               ioComms));
        }
        else
        {
          throw Exception() << "Unknown property output format '" << format << "'";
        }
      }
    }

    PropertyWriter::~PropertyWriter()
    {
      for (unsigned outputNumber = 0; outputNumber < localPropertyOutputs.size(); ++outputNumber)
      {
        delete localPropertyOutputs[outputNumber];
      }
    }

    const std::vector<LocalPropertyOutput*>& PropertyWriter::GetPropertyOutputs() const
    {
      return localPropertyOutputs;
    }

    void PropertyWriter::Write(unsigned long iterationNumber, unsigned long totalSteps) const
    {
      for (unsigned outputNumber = 0; outputNumber < localPropertyOutputs.size(); ++outputNumber)
      {
        localPropertyOutputs[outputNumber]->Write((uint64_t) iterationNumber, totalSteps);
      }
    }
  }
}
