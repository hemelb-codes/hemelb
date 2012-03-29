#include "extraction/PropertyActor.h"

namespace hemelb
{
  namespace extraction
  {
    PropertyActor::PropertyActor(const lb::SimulationState& simulationState,
                                 const std::vector<PropertyOutputFile*>& propertyOutputs,
                                 IterableDataSource& dataSource) :
        simulationState(simulationState)
    {
      propertyWriter = new PropertyWriter(dataSource, propertyOutputs);
    }

    PropertyActor::~PropertyActor()
    {
      delete propertyWriter;
    }

    void PropertyActor::SetRequiredProperties(lb::MacroscopicPropertyCache& propertyCache)
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
                propertyCache.densityCache.SetRefreshFlag();
                break;
              case OutputField::Velocity:
                propertyCache.velocityCache.SetRefreshFlag();
                break;
              case OutputField::ShearStress:
                propertyCache.shearStressCache.SetRefreshFlag();
                break;
              case OutputField::VonMisesStress:
                propertyCache.vonMisesStressCache.SetRefreshFlag();
                break;
              case OutputField::ShearRate:
                propertyCache.shearRateCache.SetRefreshFlag();
                break;
            }
          }
        }
      }
    }

    void PropertyActor::EndIteration()
    {
      propertyWriter->Write(simulationState.GetTimeStep());
    }

  }
}
