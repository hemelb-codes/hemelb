#include "extraction/PropertyActor.h"

namespace hemelb
{
  namespace extraction
  {
    PropertyActor::PropertyActor(const lb::SimulationState& simulationState,
                                 const std::vector<PropertyOutputFile>& propertyOutputs,
                                 IterableDataSource& dataSource) :
      simulationState(simulationState), propertyWriter(dataSource, propertyOutputs)
    {

    }

    void PropertyActor::EndIteration()
    {
      propertyWriter.Write(simulationState.GetTimeStep());
    }

  }
}
