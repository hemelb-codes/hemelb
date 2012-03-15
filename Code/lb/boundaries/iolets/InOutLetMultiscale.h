#ifndef HEMELB_LB_BOUNDARIES_IOLETS_INOUTLETMULTISCALE_H
#define HEMELB_LB_BOUNDARIES_IOLETS_INOUTLETMULTISCALE_H

#include "lb/boundaries/iolets/InOutLet.h"
#include "multiscale/Intercommunicand.h"
#include "multiscale/SharedValue.h"
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
                multiscale::Intercommunicand(), InOutLet(), Pressure(this), MinPressure(this), MaxPressure(this), Velocity(this,
                                                                                                                           0.0)
            {

            }
            virtual ~InOutLetMultiscale()
            {
            }
            virtual void DoIO(TiXmlElement *iParent, bool iIsLoading, configuration::SimConfig* iSimConfig);
            virtual InOutLet* Clone()
            {
              InOutLetMultiscale* copy = new InOutLetMultiscale(*this);
              return copy;
            }
            virtual void Reset(SimulationState &state)
            {
              //pass;
            }
            virtual bool IsRegistrationRequired()
            {
              return true;
            }
            LatticeDensity GetDensity(unsigned long time_step)
            {
              return mUnits->ConvertPressureToLatticeUnits(Pressure);
            }
            PhysicalPressure GetPressureMin()
            {
              return MinPressure;
            }
            PhysicalPressure GetPressureMax()
            {
              return MaxPressure;
            }

            std::string Label;

            multiscale::SharedValue<PhysicalPressure> Pressure;
            multiscale::SharedValue<PhysicalPressure> MinPressure;
            multiscale::SharedValue<PhysicalPressure> MaxPressure;
            multiscale::SharedValue<PhysicalVelocity> Velocity;

            template<class Intercommunicator> void Register(Intercommunicator &intercomms,
                                                            typename Intercommunicator::IntercommunicandTypeT &type)
            {
              intercomms.RegisterIntercommunicand(type, *this, Label);
            }
            template<class IntercommunicandType> static void DefineType(IntercommunicandType &type)
            {
              // The intercommunicators have a shared buffer which represents imaginary communication
              type.template RegisterSharedValue<PhysicalPressure>("pressure");
              type.template RegisterSharedValue<PhysicalPressure>("minPressure");
              type.template RegisterSharedValue<PhysicalPressure>("maxPressure");
              type.template RegisterSharedValue<PhysicalPressure>("velocity");
            }
        };
      }
    }
  }
}

#endif /* HEMELB_LB_BOUNDARIES_IOLETS_INOUTLETMULTISCALE_H */
