// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 
// Multiscale branch created 5/9/12
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
        virtual ~Intercommunicator()
        {

        }

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
         * Share multiscale information.
         * @param newtime time advanced to last step.
         * @returntrue if HemeLB should advance again.
         */
        virtual bool DoMultiscale(double newtime)=0;

        /***
         * Share initial multiscale information to the partners, achieving consistent initial conditions
         */
        virtual void ShareInitialConditions()=0;

        void RegisterIntercommunicand(IntercommunicandTypeT & resolver,
                                      Intercommunicand & intercommunicand,
                                      const std::string &label)
        {
          registeredObjects.insert(std::make_pair(&intercommunicand, std::make_pair(&resolver, label)));

        }

      protected:
        typedef std::map<Intercommunicand *, std::pair<IntercommunicandTypeT *, std::string> > ContentsType;
        ContentsType registeredObjects;
    };
  }
}

#endif // HEMELB_MULTISCALE_INTERCOMMUNICATOR_H
