
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "lb/iolets/InOutLetMultiscale.h"
#include "configuration/SimConfig.h"
#include "lb/iolets/BoundaryValues.h"
#include "comm/Async.h"

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

      void InOutLetMultiscale::Begin(BoundaryValues* bv) {
        //TODO: Change these operators on SharedValue.
        pressure_array[0] = pressure.GetPayload();
        pressure_array[1] = minPressure.GetPayload();
        pressure_array[2] = maxPressure.GetPayload();
      }
      void InOutLetMultiscale::Receive(BoundaryValues* bv, comm::Async::Ptr commQ) {
	// everyone receives from BC master proc
	commQ->Irecv(pressure_array, 3, bv->GetBCProcRank(), 0);
      }
      void InOutLetMultiscale::Send(BoundaryValues* bv, comm::Async::Ptr commQ) {
	if(bv->IsCurrentProcTheBCProc())
	{
	  const std::vector<int>& procList = bv->GetProcsForIolet(this); //TODO: CHECK + IMPROVE!
	  
	  for (auto dest_rank: procList)
	    commQ->Isend(pressure_array, 3, dest_rank, 0);
	}
      }
      
      void InOutLetMultiscale::CommsComplete(BoundaryValues* bv) {
	if (!bv->IsCurrentProcTheBCProc())
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
