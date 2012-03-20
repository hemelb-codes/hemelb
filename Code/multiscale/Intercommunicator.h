#ifndef HEMELB_MULTISCALE_INTERCOMMUNICATOR_H
#define HEMELB_MULTISCALE_INTERCOMMUNICATOR_H

#include "multiscale/IntercommunicandType.h"
#include "multiscale/Intercommunicand.h"
#include <map>

namespace hemelb
{
  namespace multiscale
  {
    /***
     * Abstract representation of an intercommunicator capable of communicating between processes.
     * Typical usage in a simulation step:
     *
     * if (!intercomms.ShouldAdvance()) return;
     * intercomms.GetFromMultiscale();
     * DoTimeStep();
     * intercomms.AdvanceTime(GetTime());
     * intercomms.SendToMultiscale();
     *
     * @tparam RuntimeTypeImplementation A choice of how to represent types at runtime, such as MPI_DATATYPE
     */
    template<class RuntimeTypeImplementation> class Intercommunicator
    {
      public:
        /***
         * The RuntimeTypeImplementation template parameter must have a typedef field RuntimeType, corresponding to the used for the RTT
         * And a GetType method giving the RTT value for a given type.
         * Contract example:
         * template<class T> struct ExampleTypeTraits
         * {
         * static unsigned int type;
         * };
         *
         * template<class T> unsigned int ExampleTypeTraits<T>::type = 0;
         * template<> unsigned int ExampleTypeTraits<double>::type = 1;
         * template<> unsigned int ExampleTypeTraits<int>::type = 2;
         *
         * struct ExampleRuntimeTypeImplementation
         * {
         *  typedef unsigned int RuntimeType;
         *  template<class T> static unsigned int GetType()
         *  {
         *      return ExampleTypeTraits<T>::type;
         *  }
         * };
         */
        typedef RuntimeTypeImplementation RuntimeTypeTraits;
        typedef typename RuntimeTypeImplementation::RuntimeType RuntimeType;
        /***
         * The type which should be used for specifying an intercommunicand type.
         * Example use:
         * void DefineType(IntercommunicandType &type)
         * {
         *  type.template RegisterSharedValue<PhysicalPressure>("pressure");
         *  type.template RegisterSharedValue<PhysicalVelocity>("velocity");
         }
         */
        typedef IntercommunicandType<RuntimeTypeImplementation> IntercommunicandTypeT;

        /***
         * Declare to the multiscale simulation, the value to which the preceding simulation step advanced the time.
         * @param newtime
         */
        void AdvanceTime(double newtime);
        /***
         * Ask the multiscale system, if this component should execute a simulation step, or idle waiting for brethren to catch up.
         * @return
         */
        bool ShouldAdvance(); // return false if this simulation is ahead in time, and should wait.

        /***
         * Declare an intercommunicand to the multiscale system.
         * This intercommunicand is a collection of shared values which the multiscale system should read and write
         * @param resolver A collection of labels and runtime type values describing the intercommunicand
         * @param intercommunicand The intercommunicand whose values should be read and written.
         * @param label A string label identifying the intercommunicand to the system.
         */
        void RegisterIntercommunicand(IntercommunicandTypeT & resolver,
                                      Intercommunicand & intercommunicand,
                                      const std::string &label)
        {
          registeredObjects.insert(std::make_pair(&intercommunicand, std::make_pair(&resolver, label)));
        }

        /***
         * Set the shared values in all registered intercommunicands, based on the values in sibling processes
         */
        void GetFromMultiscale();
        /***
         * Get the values from registered intercommunicands.
         */
        void SendToMultiscale();

      protected:
        typedef std::map<Intercommunicand *, std::pair<IntercommunicandTypeT *, std::string> > ContentsType;
        ContentsType registeredObjects;
    };
  }
}

#endif // HEMELB_MULTISCALE_INTERCOMMUNICATOR_H
