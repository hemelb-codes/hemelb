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

    void PropertyActor::EndIteration()
    {
      propertyWriter->Write(simulationState.GetTimeStep());
    }

  }
}
