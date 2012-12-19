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
            InOutLetMultiscale();
            /***
             * The shared values are registered through the initialiser-list.
             * pressure is initialized using other.GetPressureMax() for now, because GetPressure() is not present in the
             * parent Iolet object, which is being cloned in BoundaryValues.h. This is a somewhat ugly workaround, but it does
             * preserve the single-scale code.
             */
            InOutLetMultiscale(const InOutLetMultiscale &other);
            virtual ~InOutLetMultiscale();

            virtual void DoIO(TiXmlElement *parent, bool isLoading, configuration::SimConfig* simConfig);

            virtual InOutLet* Clone() const;
            virtual void Reset(SimulationState &state);
            virtual bool IsRegistrationRequired() const;

            LatticeDensity GetDensity(unsigned long timeStep) const;
            PhysicalPressure GetPressureMin() const;
            PhysicalPressure GetPressureMax() const;
            PhysicalVelocity GetVelocity() const;
            PhysicalPressure GetPressure() const;

            multiscale::SharedValue<PhysicalPressure> & GetPressureReference();
            PhysicalVelocity & GetVelocityReference();

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
            }

            // This should be const, and we should have a setter.
            // But the way SimConfig is set up prevents this.
            std::string & GetLabel();

            /* TODO: Be a bit smarter about this. IoletMS might not *always* require communications, but it's better now to be
             * inefficient than to end up with an inconsistent state. */
            virtual bool IsCommsRequired() const;
            virtual void SetCommsRequired(bool b);
            void DoComms(bool isIoProcess, const LatticeTime timeStep);

          protected:
            std::string label;

            bool commsRequired;
            multiscale::SharedValue<PhysicalPressure> pressure;
            multiscale::SharedValue<PhysicalPressure> minPressure;
            multiscale::SharedValue<PhysicalPressure> maxPressure;
            mutable PhysicalVelocity velocity;

        };
      }
    }
  }
}

#endif /* HEMELB_LB_BOUNDARIES_IOLETS_INOUTLETMULTISCALE_H */
