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

        namespace multiscale_constants
        {
          const PhysicalPressure HEMELB_MULTISCALE_REFERENCE_PRESSURE = 80.0;
          const PhysicalVelocity HEMELB_MULTISCALE_REFERENCE_VELOCITY = 0.0;
        }

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
                multiscale::Intercommunicand(), InOutLet(), pressure(this,
                                                                     multiscale_constants::HEMELB_MULTISCALE_REFERENCE_PRESSURE), minPressure(this,
                                                                                                                                              multiscale_constants::HEMELB_MULTISCALE_REFERENCE_PRESSURE), maxPressure(this,
                                                                                                                                                                                                                       multiscale_constants::HEMELB_MULTISCALE_REFERENCE_PRESSURE), velocity(multiscale_constants::HEMELB_MULTISCALE_REFERENCE_VELOCITY)
            {
              hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("On Creation: IoletMS, GetDensity(): velocity is %f, pressure is %f or %f",
                                                                                   GetVelocity(),
                                                                                   GetPressure(),
                                                                                   GetPressureMax());
            }
            /***
             * The shared values are registered through the initialiser-list syntactic sugar.
             * pressure is initialized using other.GetPressureMax() for now, because GetPressure() is not present in the
             * parent Iolet object, which is being cloned in BoundaryValues.h. This is a somewhat ugly workaround, but it does
             * preserve the single-scale code.
             */
            InOutLetMultiscale(const InOutLetMultiscale &other) :
                Intercommunicand(other), label(other.label), commsRequired(false), pressure(this, other.GetPressureMax()), minPressure(this,
                                                                                                                 other.GetPressureMin()), maxPressure(this,
                                                                                                                                                      other.GetPressureMax()), velocity(other.GetVelocity())
            {
              hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("On Clone: IoletMS, GetDensity(): velocity is %f, pressure is %f/%f or %f/%f",
                                                                                   GetVelocity(),
                                                                                   GetPressure(),
                                                                                   other.GetPressure(),
                                                                                   GetPressureMax(),
                                                                                   other.GetPressureMax());
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
              velocity = GetVelocity();
              //hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::OnePerCore>("IoletMS, GetDensity(): density is %f, velocity is %f, pressure is %f or %f",
              //                                                                     units->ConvertPressureToLatticeUnits(maxPressure)
              //                                                                         / Cs2,
              //                                                                     GetVelocity(),
              //                                                                     GetPressure(),
              //                                                                     GetPressureMax());
              /* TODO: Fix pressure and GetPressure values (using PressureMax() for now). */
              return units->ConvertPressureToLatticeUnits(maxPressure.GetPayload()) / Cs2;
            }
            PhysicalPressure GetPressureMin() const
            {
              return minPressure.GetPayload();
            }
            PhysicalPressure GetPressureMax() const
            {
              return maxPressure.GetPayload();
            }
            virtual PhysicalVelocity GetVelocity() const
            {
              return velocity;
            }
            PhysicalPressure GetPressure() const
            {
              return pressure.GetPayload();
            }

            multiscale::SharedValue<PhysicalPressure> & GetPressureReference()
            {
              return pressure;
            }

            //multiscale::SharedValue<PhysicalVelocity> & GetVelocityReference()
            PhysicalVelocity & GetVelocityReference()
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
              //type.template RegisterSharedValue<PhysicalVelocity>("velocity"); //leave vel out of this. Should put in IOLetVA...
            }
            // This should be const, and we should have a setter.
            // But the way SimConfig is set up prevents this.
            std::string & GetLabel()
            {
              return label;
            }

            /* TODO: Be a bit smarter about this. IoletMS might not *always* require communications, but it's better now to be
             * inefficient than to end up with an inconsistent state. */
            virtual bool IsCommsRequired() const
            {
              return commsRequired;
            }

            virtual void SetCommsRequired(bool b) {
              commsRequired = b;
            }

            void DoComms(bool isIoProcess,const LatticeTime timeStep);

          protected:
            std::string label;

            bool commsRequired;
            multiscale::SharedValue<PhysicalPressure> pressure;
            multiscale::SharedValue<PhysicalPressure> minPressure;
            multiscale::SharedValue<PhysicalPressure> maxPressure;
            mutable PhysicalVelocity velocity;
            //mutable multiscale::SharedValue<PhysicalVelocity> velocity;

        };
      }
    }
  }
}

#endif /* HEMELB_LB_BOUNDARIES_IOLETS_INOUTLETMULTISCALE_H */
