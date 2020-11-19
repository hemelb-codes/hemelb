// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_TESTS_MULTISCALE_MOCKINTERCOMMUNICATOR_H
#define HEMELB_TESTS_MULTISCALE_MOCKINTERCOMMUNICATOR_H

#include <algorithm>
#include <functional>
#include <vector>
#include <sstream>

#include "multiscale/Intercommunicator.h"
#include "net/mpi.h"

namespace hemelb
{
  namespace tests
  {
    /***
     * Example type traits structure, using the HemeLB implementation of MPI_Datatype traits.
     */

    struct MPIRuntimeType
    {
      typedef MPI_Datatype RuntimeType;

      template<class T>
      static RuntimeType GetType()
      {
	return net::MpiDataTypeTraits<T>::GetMpiDataType();
      }
    };

    /***
     * This is a very dumb example of an intercommunicator
     * It stores communicated examples in a string-keyed buffer
     * By sharing the same buffer between multiple intercommunicator interfaces, one can mock the behaviour of interprocess communication.
     * Orchestration is true if it is intent out, and false if intent in.
     */
    class MockIntercommunicator : public hemelb::multiscale::Intercommunicator<MPIRuntimeType>
    {
    public:
      MockIntercommunicator(std::map<std::string, double> & buffer, std::map<std::string, bool> &orchestration);

      void ShareInitialConditions();
      bool DoMultiscale(double new_time);

    private:
      void AdvanceTime(double new_time);
      bool ShouldAdvance();
      bool GetFromMultiscale();
      void SendToMultiscale();
      void Receive(const std::string & fieldLabel,
		   RuntimeType type,
		   const std::string objectLabel,
		   hemelb::multiscale::BaseSharedValue & value);
      void Send(const std::string & fieldLabel,
		RuntimeType type,
		const std::string objectLabel,
		hemelb::multiscale::BaseSharedValue & value);

      std::map<std::string, double> &doubleContents;
      double currentTime;
      std::map<std::string, bool> & orchestration;
    };

  }
}

#endif // HEMELB_TEST_MULTISCALE_MOCKINTERCOMMUNICATOR_H
