// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_LB_IOLETS_INOUTLETMULTISCALE_H
#define HEMELB_LB_IOLETS_INOUTLETMULTISCALE_H

#include "lb/iolets/InOutLet.h"
#include "multiscale/Intercommunicand.h"
#include "multiscale/SharedValue.h"
#include "log/Logger.h"
namespace hemelb
{
  namespace lb
  {
    namespace iolets
    {
      namespace multiscale_constants
      {
        const PhysicalPressure HEMELB_MULTISCALE_REFERENCE_PRESSURE = 0.0;
        const PhysicalVelocity HEMELB_MULTISCALE_REFERENCE_VELOCITY = PhysicalVelocity(0.0, 0.0, 0.0);
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
       * ExchangeAreaSize is the size of the area which is used to exchange information with with the outside world. It is
       * set to '1' for 1-dimensional iolets, to the surface area (in lattice sites) of the InOutLet for 2-dimensional iolets,
       * and to an even higher value should we wish to go for overlapping exchange regions.
       */
      class InOutLetMultiscale : public multiscale::Intercommunicand,
                                 public InOutLet
      {
        public:
          InOutLetMultiscale();
          InOutLetMultiscale(const InOutLetMultiscale &other);
          virtual ~InOutLetMultiscale();
          virtual InOutLet* Clone() const;
          virtual void Reset(SimulationState &state);
          virtual bool IsRegistrationRequired() const;

          virtual int GetNumberOfFieldPoints() const;
          // returns the number of field points that are exchanged with the coupled code.
          // field points will map one-to-one to HemeLB lattice sites initially; we will
          // rely on external tools to perform any interpolation tasks.
          LatticeDensity GetDensity(unsigned long timeStep) const;
          PhysicalPressure GetPressureMin() const;
          PhysicalPressure GetPressureMax() const;
          PhysicalVelocity GetVelocity() const;
          PhysicalPressure GetPressure() const;
          //std::vector<PhysicalPressure GetPressures() const;
          //TODO: update Velocity data type.
          //TODO: enable array returns in these functions (for pressure and velocity) when we have more than one field point.

          multiscale::SharedValue<PhysicalPressure> & GetPressureReference();
          multiscale::SharedValue<PhysicalVelocity> & GetVelocityReference();

          template<class Intercommunicator> void Register(
                                                          Intercommunicator &intercomms,
                                                          typename Intercommunicator::IntercommunicandTypeT &type)
          {
            intercomms.RegisterIntercommunicand(type, *this, label);
          }

          template<class IntercommunicandType> static void DefineType(IntercommunicandType &type)
          {
            // The intercommunicators have a shared buffer which represents imaginary communication
            type.template RegisterSharedValue<PhysicalPressure> ("pressure");
            type.template RegisterSharedValue<PhysicalPressure> ("minPressure");
            type.template RegisterSharedValue<PhysicalPressure> ("maxPressure");
            //type.template RegisterSharedValue<PhysicalPressure>("velocity");
          }
          // This should be const, and we should have a setter.
          // But the way SimConfig is set up prevents this.
          std::string & GetLabel();

          virtual bool IsCommsRequired() const;
          virtual void SetCommsRequired(bool b);
          void DoComms(bool isIoProcess, const LatticeTime timeStep);

        private:
          std::string label;

          bool commsRequired;
          int numberOfFieldPoints;
          multiscale::SharedValue<PhysicalPressure> pressure;
          multiscale::SharedValue<PhysicalPressure> minPressure;
          multiscale::SharedValue<PhysicalPressure> maxPressure;
          mutable multiscale::SharedValue<PhysicalVelocity> velocity;
      };
    }
  }
}

#endif /* HEMELB_LB_IOLETS_INOUTLETMULTISCALE_H */
