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
      class BoundaryComms;
      //forward decl;
      namespace iolets
      {

        class InOutLet
        {
          public:
            InOutLet();
            virtual ~InOutLet();

            /***
             * Read the TinyXML structure and set up the iolet, or write it to a TinyXML structure.
             * @param iParent Parent XML element
             * @param iIsLoading Read or write?
             * @param simConfig The config object being read
             */
            virtual void DoIO(TiXmlElement *iParent, bool iIsLoading, configuration::SimConfig* simConfig) = 0;

            /***
             * Copy the InOutLet.
             * @return Pointer to new IOLet.
             */
            virtual InOutLet* Clone() = 0;

            /***
             * Carry out any communication needed for the IOLet.
             * @return true if any comms were done.
             */
            virtual bool GetIsCommsRequired()
            {
              return false;
            }

            void SetComms(BoundaryComms * acomms)
            {
              comms = acomms;
            }
            BoundaryComms * GetComms()
            {
              return comms;
            }
            virtual void DoComms(bool is_io_proc)
            {
              // pass
            }
            /***
             * Set up the Iolet.
             * @param units a UnitConverter instance.
             */
            void Initialise(const util::UnitConverter* units);

            /***
             * Get the minimum density, in lattice units
             * @return minimum density, in lattice units
             */
            LatticePressure GetDensityMin();

            /***
             * Get the maximum density, in lattice units
             * @return maximum density, in lattice units
             */
            LatticePressure GetDensityMax();

            /***
             * Get the minimum pressure, in physical units
             * @return
             */
            virtual PhysicalPressure GetPressureMin()=0;

            /***
             * Get the maximum pressure, in physical units
             * @return
             */
            virtual PhysicalPressure GetPressureMax()=0;
            virtual LatticeDensity GetDensity(LatticeTime time_step)=0;
            virtual void Reset(SimulationState &state)=0;
            util::Vector3D<float> Position; //! !!!!! Public data member!
            util::Vector3D<float> Normal; //! !!!!!! Public data member!
          protected:
            const util::UnitConverter* mUnits;
            BoundaryComms * comms;
        };
      }
    }
  }
}

#endif /* HEMELB_LB_BOUNDARIES_IOLETS_INOUTLET_H */
