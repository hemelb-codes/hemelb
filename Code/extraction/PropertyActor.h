#ifndef HEMELB_EXTRACTION_PROPERTYACTOR_H
#define HEMELB_EXTRACTION_PROPERTYACTOR_H

#include "extraction/PropertyWriter.h"
#include "lb/SimulationState.h"
#include "net/IteratedAction.h"

namespace hemelb
{
  namespace extraction
  {
    class PropertyActor : public net::IteratedAction
    {
      public:
        /**
         * Constructor, gets the class ready for reading.
         * @param simulationState
         * @param propertyOutputs
         * @param dataSource
         * @return
         */
        PropertyActor(const lb::SimulationState& simulationState,
                      const std::vector<PropertyOutputFile>& propertyOutputs,
                      IterableDataSource& dataSource);

        /**
         * Override the iterated actor end of iteration method to perform writing.
         */
        void EndIteration();

      private:
        const lb::SimulationState& simulationState;
        PropertyWriter propertyWriter;
    };
  }
}

#endif /* HEMELB_EXTRACTION_PROPERTYACTOR_H */
