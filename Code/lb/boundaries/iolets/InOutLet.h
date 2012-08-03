// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_LB_BOUNDARIES_IOLETS_INOUTLET_H
#define HEMELB_LB_BOUNDARIES_IOLETS_INOUTLET_H

#include "util/Vector3D.h"
#include "util/UnitConverter.h"
#include "tinyxml.h"

namespace hemelb
{
  namespace configuration
  {
    class SimConfig;
  }

  namespace lb
  {
    namespace boundaries
    {

      //forward declare boundary comms class
      class BoundaryComms;
      namespace iolets
      {
        /**
         * Base Iolet class
         * Contains information configured from the xml config file, and calculates a density near itself for use in LB calculation
         * Provides maximum and minimum range of densities/pressures for use by steering.
         */
        class InOutLet
        {
          public:
            InOutLet() :
                comms(NULL)
            {
            }
            virtual ~InOutLet()
            {
            }

            /***
             * Read the TinyXML structure and set up the iolet, or write it to a TinyXML structure.
             * @param iParent Parent XML element
             * @param iIsLoading Read or write?
             * @param simConfig The config object being read
             */
            virtual void DoIO(TiXmlElement *parent, bool isLoading, configuration::SimConfig* simConfig) = 0;

            /***
             * Copy the InOutLet.
             * @return Pointer to new IOLet.
             */
            virtual InOutLet* Clone() const = 0;

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
            virtual void DoComms(bool isIoProcess,const LatticeTime timeStep)
            {
              // pass
            }
            /***
             * Set up the Iolet.
             * @param units a UnitConverter instance.
             */
            void Initialise(const util::UnitConverter* unitConverter)
            {
              units = unitConverter;
            }

            /***
             * Get the minimum density, in lattice units
             * @return minimum density, in lattice units
             */
            LatticePressure GetDensityMin() const
            {
              return units->ConvertPressureToLatticeUnits(GetPressureMin()) / Cs2;
            }

            /***
             * Get the maximum density, in lattice units
             * @return maximum density, in lattice units
             */
            LatticePressure GetDensityMax() const
            {
              return units->ConvertPressureToLatticeUnits(GetPressureMax()) / Cs2;
            }

            /***
             * Get the minimum pressure, in physical units
             * @return
             */
            virtual PhysicalPressure GetPressureMin() const =0;

            /***
             * Get the maximum pressure, in physical units
             * @return
             */
            virtual PhysicalPressure GetPressureMax() const =0;
            virtual LatticeDensity GetDensity(LatticeTime time_step) const =0;
            virtual void Reset(SimulationState &state)=0;
            // TODO I do not like returning non-const references, this method should be const and we should have a setter.
            // but, the way the IO code in SimConfig is currently set up prevents this for now.
            util::Vector3D<float> &GetPosition()
            {
              return position;
            }
            // TODO I do not like returning non-const references, this method should be const and we should have a setter.
            // but, the way the IO code in SimConfig is currently set up prevents this for now.
            util::Vector3D<float> &GetNormal()
            {
              return normal;
            }
          protected:
            util::Vector3D<float> position;
            util::Vector3D<float> normal;
            const util::UnitConverter* units;
            BoundaryComms * comms;
        };
      }
    }
  }
}

#endif /* HEMELB_LB_BOUNDARIES_IOLETS_INOUTLET_H */
