// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#include "lb/iolets/InOutLetMultiscale.h"
#include "configuration/SimConfig.h"
#include "net/IOCommunicator.h"
#include "lb/iolets/BoundaryComms.h"
#include "lb/iolets/BoundaryValues.h"

namespace hemelb
{
  namespace lb
  {
    namespace iolets
    {

      InOutLetMultiscale::InOutLetMultiscale() :
        multiscale::Intercommunicand(), InOutLet(),
            units(NULL),
            numberOfFieldPoints(1),
            pressure(this, multiscale_constants::HEMELB_MULTISCALE_REFERENCE_PRESSURE),
            minPressure(this, multiscale_constants::HEMELB_MULTISCALE_REFERENCE_PRESSURE),
            maxPressure(this, multiscale_constants::HEMELB_MULTISCALE_REFERENCE_PRESSURE),
            velocity(this, multiscale_constants::HEMELB_MULTISCALE_REFERENCE_VELOCITY)
      {
      }
      /***
       * The shared values are registered through the initialiser-list syntactic sugar.
       * pressure is initialized using other.GetPressureMax() for now, because GetPressure() is not present in the
       * parent Iolet object, which is being cloned in BoundaryValues.h. This is a somewhat ugly workaround, but it does
       * preserve the single-scale code.
       */
      InOutLetMultiscale::InOutLetMultiscale(const InOutLetMultiscale &other) :
        Intercommunicand(other), label(other.label), units(other.units), commsRequired(false),
            pressure(this, other.maxPressure.GetPayload()), minPressure(this, other.minPressure.GetPayload()),
            maxPressure(this, other.maxPressure.GetPayload()), velocity(this, other.GetVelocity())
      {
      }

      void InOutLetMultiscale::Initialise(const util::UnitConverter* unitConverter)
      {
        units = unitConverter;
      }

      InOutLetMultiscale::~InOutLetMultiscale()
      {
      }

      InOutLet* InOutLetMultiscale::Clone() const
      {
        InOutLetMultiscale* copy = new InOutLetMultiscale(*this);
        return copy;
      }
      void InOutLetMultiscale::Reset(SimulationState &state)
      {
        //pass;
      }
      bool InOutLetMultiscale::IsRegistrationRequired() const
      {
        return true;
      }
      int InOutLetMultiscale::GetNumberOfFieldPoints() const
      {
        return numberOfFieldPoints;
      }
      LatticeDensity InOutLetMultiscale::GetDensity(unsigned long timeStep) const
      {
        /* TODO: Fix pressure and GetPressure values (using PressureMax() for now). */
        return units->ConvertPressureToLatticeUnits(maxPressure.GetPayload()) / Cs2;
      }
      LatticeDensity InOutLetMultiscale::GetDensityMin() const
      {
        return units->ConvertPressureToLatticeUnits(minPressure.GetPayload()) / Cs2;
      }
      LatticeDensity InOutLetMultiscale::GetDensityMax() const
      {
        return units->ConvertPressureToLatticeUnits(maxPressure.GetPayload()) / Cs2;
      }
      PhysicalVelocity InOutLetMultiscale::GetVelocity() const
      {
        return velocity;
      }
      PhysicalPressure InOutLetMultiscale::GetPressure() const
      {
        return pressure.GetPayload();
      }

      multiscale::SharedValue<PhysicalPressure> & InOutLetMultiscale::GetPressureReference()
      {
        return pressure;
      }

      multiscale::SharedValue<PhysicalVelocity> & InOutLetMultiscale::GetVelocityReference()
      {
        return velocity;
      }

      // This should be const, and we should have a setter.
      // But the way SimConfig is set up prevents this.
      std::string & InOutLetMultiscale::GetLabel()
      {
        return label;
      }

      /* TODO: Be a bit smarter about this. IoletMS might not *always* require communications, but it's better now to be
       * inefficient than to end up with an inconsistent state. */
      bool InOutLetMultiscale::IsCommsRequired() const
      {
        return commsRequired;
      }

      void InOutLetMultiscale::SetCommsRequired(bool b)
      {
        commsRequired = b;
      }

      /* Distribution of internal pressure values */
      void InOutLetMultiscale::DoComms(bool isIoProc, LatticeTimeStep time_step)
      {
        hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("DoComms in IoletMultiscale triggered: %s",
                                                                              isIoProc
                                                                                ? "true"
                                                                                : "false");
        double pressure_array[3];
        //TODO: Change these operators on SharedValue.
        pressure_array[0] = pressure.GetPayload();
        pressure_array[1] = minPressure.GetPayload();
        pressure_array[2] = maxPressure.GetPayload();

        net::Net commsNet;

        const std::vector<int>& procList = comms->GetListOfProcs(); //TODO: CHECK + IMPROVE!

        // If this proc is to do IO, send the pressure array list to all cores that require it.
        if (isIoProc && procList[0] != BoundaryValues::GetBCProcRank())
        {
          for (std::vector<int>::const_iterator it = procList.begin(); it != procList.end(); it++)
          {
            commsNet.RequestSend(pressure_array, 3, *it);
          }
        }
        // Otherwise, receive the pressure array list from the core.
        else if (procList[0] != BoundaryValues::GetBCProcRank())
        {
          commsNet.RequestReceive(pressure_array, 3, BoundaryValues::GetBCProcRank());
        }

        // Perform the send / receive.
        commsNet.Dispatch();

        if (!isIoProc)
        {
          pressure.SetPayload(static_cast<PhysicalPressure> (pressure_array[0]));
          minPressure.SetPayload(static_cast<PhysicalPressure> (pressure_array[1]));
          maxPressure.SetPayload(static_cast<PhysicalPressure> (pressure_array[2]));
          hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("Received: %f %f %f",
                                                                                pressure.GetPayload(),
                                                                                minPressure.GetPayload(),
                                                                                maxPressure.GetPayload());
        }
      }
    }
  }
}
