// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

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

        /***
         * An inlet/outlet whose density is obtained from the multiscale intercommunicator
         * We envisage communication of velocity information outwards to other processes
         * the velocity SharedValue is a place-holder for this.
         * We set 0.1 as an arbitrary value for the output velocity for now.
         * We do not yet have an understanding of the necessary physics, or the necessary computational infrastructure,
         * for building boundary conditions which have access to their nearby velocity field.
         * The min and max pressure SharedValues are placeholders for information needed by the steering and visualisation code.
         * The steering and visualisation code requires minimum and maximum pressure values.
         */
        class InOutLetMultiscale : public multiscale::Intercommunicand,
                                   public InOutLet
        {
          public:
            InOutLetMultiscale() :
                multiscale::Intercommunicand(), InOutLet(), pressure(this), minPressure(this), maxPressure(this), velocity(this)
            {

            }
            /***
             * The shared values are registered through the initialiser-list syntactic sugar.
             */
            InOutLetMultiscale(const InOutLetMultiscale &other) :
                Intercommunicand(other), label(other.label), pressure(this, other.GetPressure()), minPressure(this,
                                                                                                              other.GetPressureMin()), maxPressure(this,
                                                                                                                                                   other.GetPressureMax()), velocity(this,
                                                                                                                                                                                     other.GetVelocity())
            {

            }
            virtual ~InOutLetMultiscale()
            {
            }
            virtual void DoIO(TiXmlElement *parent, bool isLoading, configuration::SimConfig* simConfig);

            virtual InOutLet* Clone() const
            {
              InOutLetMultiscale* copy = new InOutLetMultiscale(*this);
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
              return units->ConvertPressureToLatticeUnits(pressure) / Cs2;
            }
            PhysicalPressure GetPressureMin() const
            {
              return minPressure;
            }
            PhysicalPressure GetPressureMax() const
            {
              return maxPressure;
            }
            PhysicalVelocity_deprecated GetVelocity() const
            {
              return velocity;
            }
            PhysicalPressure GetPressure() const
            {
              return pressure;
            }

            multiscale::SharedValue<PhysicalPressure> & GetPressureReference()
            {
              return pressure;
            }

            multiscale::SharedValue<PhysicalVelocity_deprecated> & GetVelocityReference()
            {
              return velocity;
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
            mutable multiscale::SharedValue<PhysicalVelocity_deprecated> velocity;
        };
      }
    }
  }
}

#endif /* HEMELB_LB_BOUNDARIES_IOLETS_INOUTLETMULTISCALE_H */
