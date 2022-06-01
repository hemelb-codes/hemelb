// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "extraction/PropertyActor.h"

namespace hemelb
{
  namespace extraction
  {
    PropertyActor::PropertyActor(const lb::SimulationState& simulationState,
                                 const std::vector<PropertyOutputFile>& propertyOutputs,
                                 IterableDataSource& dataSource, reporting::Timers& timers,
                                 const net::IOCommunicator& ioComms) :
        simulationState(simulationState), timers(timers)
    {
      propertyWriter = new PropertyWriter(dataSource, propertyOutputs, ioComms);
    }

    PropertyActor::~PropertyActor()
    {
      delete propertyWriter;
    }

    void PropertyActor::SetRequiredProperties(lb::MacroscopicPropertyCache& propertyCache)
    {
      const std::vector<LocalPropertyOutput*>& propertyOutputs =
          propertyWriter->GetPropertyOutputs();

      // Iterate over each property output spec.
      for (unsigned output = 0; output < propertyOutputs.size(); ++output)
      {
        const LocalPropertyOutput* propertyOutput = propertyOutputs[output];

        // Only consider the ones that are being written this iteration.
        if (propertyOutput->ShouldWrite(simulationState.GetTimeStep()))
        {
          auto& outputFile = propertyOutput->GetOutputSpec();

          // Iterate over each field.
	  for (auto&& fieldSpec: outputFile.fields)
          {
            // Set the cache to calculate each required field.
	    source::visit(
	      fieldSpec.src,
	      [&](source::Pressure) {
		propertyCache.densityCache.SetRefreshFlag();
	      },
	      [&](source::Velocity) {
		propertyCache.velocityCache.SetRefreshFlag();
	      },
	      [&](source::ShearStress) {
		propertyCache.wallShearStressMagnitudeCache.SetRefreshFlag();
	      },
	      [&](source::VonMisesStress) {
		propertyCache.vonMisesStressCache.SetRefreshFlag();
	      },
	      [&](source::ShearRate) {
		propertyCache.shearRateCache.SetRefreshFlag();
	      },
	      [&](source::StressTensor) {
		propertyCache.stressTensorCache.SetRefreshFlag();
	      },
	      [&](source::Traction) {
		propertyCache.tractionCache.SetRefreshFlag();
	      },
	      [&](source::TangentialProjectionTraction) {
		propertyCache.tangentialProjectionTractionCache.SetRefreshFlag();
	      },
	      [](source::Distributions) {
		// We don't actually have to cache anything to get the distribution.
	      },
	      [](source::MpiRank) {
		// We don't actually have to cache anything to get the rank.
	      }
	    );

          }
        }
      }
    }

    void PropertyActor::EndIteration()
    {
      timers[reporting::Timers::extractionWriting].Start();
      propertyWriter->Write(simulationState.GetTimeStep());
      timers[reporting::Timers::extractionWriting].Stop();
    }

  }
}
