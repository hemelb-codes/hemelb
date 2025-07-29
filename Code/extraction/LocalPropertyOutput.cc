// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "extraction/LocalPropertyOutput.h"
#include "net/IOCommunicator.h"
#include "util/variant.h"


namespace hemelb
{
  namespace extraction
  {
    // basic constructor for LocalPropertyOutput
    LocalPropertyOutput::LocalPropertyOutput(IterableDataSource& dataSource,
                                             const PropertyOutputFile& outputSpec_,
                                             const net::IOCommunicator& ioComms) :
        comms(ioComms), dataSource(dataSource), outputSpec(outputSpec_)
    {
		// Calculate local and global write points
		local_site_count = CountWrittenSitesOnRank();
		global_site_count = comms.AllReduce(local_site_count, MPI_SUM);
    }
    // write check function
    bool LocalPropertyOutput::ShouldWrite(unsigned long timestepNumber) const
    {
    	return ( (timestepNumber % outputSpec.frequency) == 0);
    }
    
    // output spec getter
    const PropertyOutputFile& LocalPropertyOutput::GetOutputSpec() const
    {
		return outputSpec;
    }

    // count the number of sites written on this rank
    uint64_t LocalPropertyOutput::CountWrittenSitesOnRank() {
      uint64_t n = 0;
      dataSource.Reset();
      while (dataSource.ReadNext())
      {
        if (outputSpec.geometry->Include(dataSource, dataSource.GetPosition()))
        {
        	++n;
        }
      }
      return n;
    }

    // get the length of the field based on the source type
	unsigned LocalPropertyOutput::GetFieldLength(source::Type src) const
    {
      return overload_visit(src,
	[](source::Pressure) {
	  return 1U;
	},
	[](source::Velocity) {
	  return 3U;
	},
	[](source::ShearStress) {
	  return 1U;
	},
	[](source::VonMisesStress) {
	  return 1U;
	},
	[](source::ShearRate) {
	  return 1U;
	},
	[](source::StressTensor) {
	  return 6U;
	},
	[](source::Traction) {
	  return 3U;
	},
	[](source::TangentialProjectionTraction) {
	  return 3U;
	},
	[&](source::Distributions) {
	  return dataSource.GetNumVectors();
	},
	[](source::MpiRank) {
	  return 1U;
	}
	  );
	}
  }
}