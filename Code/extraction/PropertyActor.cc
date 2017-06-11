
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
                                 const std::vector<PropertyOutputFile*>& propertyOutputs,
                                 std::map<OutputField::FieldType, IterableDataSource**>& dataSources,
                                 reporting::Timers& timers,
                                 const net::IOCommunicator& ioComms) :
        simulationState(simulationState), timers(timers)
    {
      propertyWriter = new PropertyWriter(dataSources, propertyOutputs, ioComms);
    }

    PropertyActor::~PropertyActor()
    {
      delete propertyWriter;
    }

    void PropertyActor::SetRequiredProperties(std::map<hemelb::extraction::OutputField::FieldType, hemelb::extraction::IterableDataSource**> dataSourceMap) // (PropertyCache& propertyCache)
    {
      const std::vector<LocalPropertyOutput*>& propertyOutputs = propertyWriter->GetPropertyOutputs();

      // Iterate over each property output spec.
      for (unsigned output = 0; output < propertyOutputs.size(); ++output)
      {
        const LocalPropertyOutput* propertyOutput = propertyOutputs[output];

        // Only consider the ones that are being written this iteration.
        if (propertyOutput->ShouldWrite(simulationState.GetTimeStep()))
        {
          const PropertyOutputFile* outputFile = propertyOutput->GetOutputSpec();

          // Iterate over each field.
          for (unsigned outputField = 0; outputField < outputFile->fields.size(); ++outputField)
          {
            // Set the cache to calculate each required field.
            switch (outputFile->fields[outputField].type)
            {
              case (OutputField::Pressure):
		(*dataSourceMap[outputFile->fields[outputField].type])->GetPropertyCache().densityCache.SetRefreshFlag();
                break;
	      case OutputField::TracerConcentration:
              case OutputField::Velocity:
		(*dataSourceMap[outputFile->fields[outputField].type])->GetPropertyCache().velocityCache.SetRefreshFlag();
		break;
              case OutputField::ShearStress:
                (*dataSourceMap[outputFile->fields[outputField].type])->GetPropertyCache().wallShearStressMagnitudeCache.SetRefreshFlag();
                break;
              case OutputField::VonMisesStress:
                (*dataSourceMap[outputFile->fields[outputField].type])->GetPropertyCache().vonMisesStressCache.SetRefreshFlag();
                break;
              case OutputField::ShearRate:
                (*dataSourceMap[outputFile->fields[outputField].type])->GetPropertyCache().shearRateCache.SetRefreshFlag();
                break;
              case OutputField::StressTensor:
                (*dataSourceMap[outputFile->fields[outputField].type])->GetPropertyCache().stressTensorCache.SetRefreshFlag();
                break;
              case OutputField::Traction:
                (*dataSourceMap[outputFile->fields[outputField].type])->GetPropertyCache().tractionCache.SetRefreshFlag();
                break;
              case OutputField::TangentialProjectionTraction:
                (*dataSourceMap[outputFile->fields[outputField].type])->GetPropertyCache().tangentialProjectionTractionCache.SetRefreshFlag();
                break;
              case OutputField::MpiRank:
                // We don't actually have to cache anything to get the rank.
                break;
              default:
                // This assert should never trip. It only occurs when someone adds a new field to OutputField
                // and forgets adding a new case to the switch
                assert(false);
            }
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
