
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_STENTS_STENT_H
#define HEMELB_LB_STENTS_STENT_H

#include "util/Vector3D.h"
#include "util/UnitConverter.h"
#include "lb/SimulationState.h"
#include "lb/iolets/BoundaryComms.h"
#include "lb/iolets/BoundaryCommunicator.h"

namespace hemelb
{
  namespace lb
  {
    namespace stents
    {

      //forward declare boundary comms class
      class BoundaryComms;
      class BoundaryCommunicator;

      /**
       * Base class for extra data needed by LB BC implementations.
       * Makes "Iolet coordinates" available.
       * These are coordinates in a frame aligned with the iolet plane.
       * iolet(0, 0, 0) corresponds to the iolet's position in lattice coordinates
       * x & y are arbitrary in plane components,
       * z is in the direction of the iolet's normal.
       */
      class Stent;

      /**
       * Base Iolet class
       * Contains information configured from the xml config file, and calculates a density near itself for use in LB calculation
       * Provides maximum and minimum range of densities/pressures for use by steering.
       */
      class Stent
      {
        public:
          Stent() :
            comms(NULL)
          {
          }
          virtual ~Stent()
          {
          }

          /***
           * Copy the InOutLet.
           * @return Pointer to new IOLet.
           */
          virtual Stent* Clone() const = 0;

          /***
           * This is a castable? virtual method, which is perhaps an anti-pattern
           * We should potentially use dynamic cast checks instead.
           * @return true if any comms were done.
           */
          virtual bool IsCommsRequired() const
          {
            return false;
          }

          /***
           * This is a castable? virtual method, which is perhaps an anti-pattern
           * We should potentially use dynamic cast checks instead.
           * @return true if any comms were done.
           */
          virtual bool IsRegistrationRequired() const
          {
            return false;
          }
          void SetComms(BoundaryComms * boundaryComms)
          {
            comms = boundaryComms;
          }
          BoundaryComms * GetComms() const
          {
            return comms;
          }
          /***
           * Carry out communication necessary
           * @param isIoProcess Is the process the master process?
           */
          virtual void DoComms(const BoundaryCommunicator& bcComms, const LatticeTimeStep timeStep);

          /***
           * Set up the Iolet.
           * @param units a UnitConverter instance.
           */
          virtual void Initialise(const util::UnitConverter* unitConverter)
          {
          }

          /***
           * Get the minimum density, in lattice units
           * @return minimum density, in lattice units
           */
          virtual LatticeDensity GetDensityMin() const = 0;

          /***
           * Get the maximum density, in lattice units
           * @return maximum density, in lattice units
           */
          virtual LatticeDensity GetDensityMax() const = 0;

          /***
           * Get the minimum pressure, in lattice units
           * @return
           */

          /***
           * Get the maximum pressure, in lattice units
           * @return
           */

          /// @todo: #632 This method must be moved to InOutletPressure
          virtual LatticeDensity GetDensity(LatticeTimeStep time_step) const = 0;

          /// @todo: #632 Is this method ever implemented not empty?
          virtual void Reset(SimulationState& state) = 0;

          /**
           * Set the normal of the InOutlet
           * @param newNormal
           */

          /**
           * Set the minimum density throughout the simulation.
           * @param minSimDensity
           */
          void SetMinimumSimulationDensity(LatticeDensity minSimDensity)
          {
            minimumSimulationDensity = minSimDensity;
          }

        protected:
          LatticeDensity minimumSimulationDensity;
          BoundaryComms* comms;
      };

    }
  }
}

#endif /* HEMELB_LB_STENTS_STENT_H */
