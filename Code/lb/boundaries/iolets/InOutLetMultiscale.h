#ifndef HEMELB_LB_BOUNDARIES_IOLETS_INOUTLETMULTISCALE_H
#define HEMELB_LB_BOUNDARIES_IOLETS_INOUTLETMULTISCALE_H

#include "lb/boundaries/iolets/InOutLet.h"
#include "multiscale/Intercommunicand.h"
#include "multiscale/SharedValue.h"
#include "log/Logger.h"
namespace hemelb
{
  namespace lb
  {
    namespace boundaries
    {
      namespace iolets
      {

        /*
         * Template values are chosen to be tUpdatePeriod = 1 and tComms = false
         * The multiscale pressure trace can be easily calculated locally by any proc and
         * there is no need to store values for any time step beyond the current one
         */
        class InOutLetMultiscale : public multiscale::Intercommunicand,
                                   public InOutLet
        {
          public:
            InOutLetMultiscale() :
                multiscale::Intercommunicand(), InOutLet(), pressure(this), minPressure(this), maxPressure(this), velocity(this,
                                                                                                                           0.0)
            {

            }
            virtual ~InOutLetMultiscale()
            {
            }
            virtual void DoIO(TiXmlElement *parent, bool isLoading, configuration::SimConfig* simConfig);
            virtual InOutLet* Clone() const
            {
              InOutLetMultiscale* copy = new InOutLetMultiscale(*this);
              // copy constructor must copy pointers to the new locations of the shared data into the copy's reference buffer.
              copy->Values().clear();
              copy->RegisterSharedValue(&copy->pressure);
              copy->RegisterSharedValue(&copy->minPressure);
              copy->RegisterSharedValue(&copy->maxPressure);
              copy->RegisterSharedValue(&copy->velocity);
              return copy;
            }
            virtual void Reset(SimulationState &state)
            {
              //pass;
            }
            virtual bool IsRegistrationRequired() const
            {
              return true;
            }
            LatticeDensity GetDensity(unsigned long timeStep) const
            {
              return units->ConvertPressureToLatticeUnits(pressure);

            }
            PhysicalPressure GetPressureMin() const
            {
              return minPressure;
            }
            PhysicalPressure GetPressureMax() const
            {
              return maxPressure;
            }

            template<class Intercommunicator> void Register(Intercommunicator &intercomms,
                                                            typename Intercommunicator::IntercommunicandTypeT &type)
            {
              intercomms.RegisterIntercommunicand(type, *this, label);
            }
            template<class IntercommunicandType> static void DefineType(IntercommunicandType &type)
            {
              // The intercommunicators have a shared buffer which represents imaginary communication
              type.template RegisterSharedValue<PhysicalPressure>("pressure");
              type.template RegisterSharedValue<PhysicalPressure>("minPressure");
              type.template RegisterSharedValue<PhysicalPressure>("maxPressure");
              type.template RegisterSharedValue<PhysicalPressure>("velocity");
            }
            // This should be const, and we should have a setter.
            // But the way SimConfig is set up prevents this.
            std::string & GetLabel()
            {
              return label;
            }
          private:
            std::string label;

            multiscale::SharedValue<PhysicalPressure> pressure;
            multiscale::SharedValue<PhysicalPressure> minPressure;
            multiscale::SharedValue<PhysicalPressure> maxPressure;
            multiscale::SharedValue<PhysicalVelocity> velocity;
        };
      }
    }
  }
}

#endif /* HEMELB_LB_BOUNDARIES_IOLETS_INOUTLETMULTISCALE_H */
