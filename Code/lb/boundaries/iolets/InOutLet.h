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

            // Should be called before simulation starts running (including after a reset)
            // Resizes densityCycle and calls CalculateCycle
            /***
             * Size the densityCycle according to the update period for this iolet,
             * and fill it in with an initial calculation.
             * @param densityCycle
             * @param iState
             */
            virtual void InitialiseCycle(std::vector<distribn_t> &densityCycle, const SimulationState *iState) = 0;
            /***
             * Determine if it is time to update the current set of stored densities (densityCycle) for this and a few subsequent steps,
             * and if so, do so by calling CalculateCycle
             * @param densityCycle
             * @param iState
             */
            virtual void UpdateCycle(std::vector<distribn_t> &densityCycle, const SimulationState *iState) = 0;
            /***
             * Fill in a vector of values of the density for this IOLet.
             * The Iolet will not subsequently calculate it's density values, but will read them from this buffer.
             * This array is NOT typically one value for each point in the pulsatile cycle.
             * The length of the densityCycle argument is an arbitrary choice of how often to calculate the densities for this and some subsequent steps.
             * This length is usually one, indicating that the calculation is done every time step,
             * or zero, indicating that the calculation is done once for the whole simulation.
             * @param densityCycle An array of densities for this iolet
             * @param iState Simulation state for the iolet
             */
            virtual void CalculateCycle(std::vector<distribn_t> &densityCycle, const SimulationState *iState) = 0;
            /***
             * Carry out any communication needed for the IOLet.
             * @return true if any comms were done.
             */
            virtual bool DoComms() = 0;

            /***
             * Reset the maximum and minimum values of the density in lattice units
             * Need to call if the simulation state has changed its voxel size or time step.
             */
            virtual void ResetValues();

            /***
             * Reset the maximum and minimum values of the density in lattice units
             * Need to call if the simulation state has changed its voxel size or time step.
             */
            void ResetCommonLatticeValues();

            /***
             * Set up the Iolet.
             * @param units a UnitConverter instance.
             */
            void Initialise(const util::UnitConverter* units);

            /***
             * Get the minimum density, in lattice units
             * @return minimum density, in lattice units
             */
            distribn_t GetDensityMin();

            /***
             * Get the maximum density, in lattice units
             * @return maximum density, in lattice units
             */
            distribn_t GetDensityMax();

            /***
             * !!!!! Public data member
             * The density at the Iolet in lattice units
             */
            distribn_t density;

            /***
             * !!!!! Public data member
             * The update period of the Iolet.
             * This is not the same as the cycle period.
             * It is the value of how often to recalculate the densities
             * !!!!! Unused, is now a template parameter in InOutLetCycle
             */
            unsigned long UpdatePeriod;

            /***
             * !!!!! Public data member
             */
            double PressureMinPhysical;
            /***
             * !!!!! Public data member
             */
            double PressureMaxPhysical;

            util::Vector3D<float> Position; //! !!!!! Public data member!
            util::Vector3D<float> Normal; //! !!!!!! Public data member!

          protected:
            const util::UnitConverter* mUnits;

          private:
            distribn_t DensityMinLattice; //! Minimum density, lattice units
            distribn_t DensityMaxLattice; //! Maximum density, lattice units

        };

      }
    }
  }
}

#endif /* HEMELB_LB_BOUNDARIES_IOLETS_INOUTLET_H */
