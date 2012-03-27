#ifndef HEMELB_UNITTESTS_MULTISCALE_MOCKINTERCOMMUNICATOR_H
#define HEMELB_UNITTESTS_MULTISCALE_MOCKINTERCOMMUNICATOR_H
#include "multiscale/Intercommunicator.h"
#include "mpiInclude.h"
#include "tinyxml.h"

#include <algorithm>
#include <functional>
#include <vector>
#include <sstream>
namespace hemelb
{
  namespace unittests
  {
    namespace multiscale
    {
      /***
       * Example type traits structure, using the HemeLB implementation of MPI_TYPE traits.
       */

      struct MPIRuntimeType
      {
          typedef MPI_Datatype RuntimeType;
          template<class T> static RuntimeType GetType()
          {
            return MpiDataTypeTraits<T>::GetMpiDataType();
          }
      };

      /***
       * This is a very dumb example of an intercommunicator
       * It stores communicated examples in a string-keyed buffer
       * By sharing the same buffer between multiple intercommunicator interfaces, one can mock the behaviour of interprocess communication.
       */
      class MockIntercommunicator : public hemelb::multiscale::Intercommunicator<MPIRuntimeType>
      {
        public:
          MockIntercommunicator(std::map<std::string, double> & buffer) :
              doubleContents(buffer), currentTime(0)
          {

          }

          void AdvanceTime(double new_time)
          {
            currentTime = new_time;
            doubleContents["shared_time"] = new_time;
          }
          bool ShouldAdvance()
          {
            return doubleContents["shared_time"] >= currentTime;
          }

          bool GetFromMultiscale()
          {
            for (ContentsType::iterator intercommunicandData = registeredObjects.begin();
                intercommunicandData != registeredObjects.end(); intercommunicandData++)
            {
              hemelb::multiscale::Intercommunicand &sharedObject = *intercommunicandData->first;
              std::string &label = intercommunicandData->second.second;
              IntercommunicandTypeT &resolver = *intercommunicandData->second.first;

              for (unsigned int sharedFieldIndex = 0; sharedFieldIndex <= sharedObject.Values().size();
                  sharedFieldIndex++)
              {
                Receive(resolver.Fields()[sharedFieldIndex].first,
                        resolver.Fields()[sharedFieldIndex].second,
                        label,
                        *sharedObject.Values()[sharedFieldIndex]);
              }
            }
            return true;
          }
          void SendToMultiscale()
          {
            for (ContentsType::iterator intercommunicandData = registeredObjects.begin();
                intercommunicandData != registeredObjects.end(); intercommunicandData++)
            {
              hemelb::multiscale::Intercommunicand &sharedObject = *intercommunicandData->first;
              std::string &label = intercommunicandData->second.second;
              IntercommunicandTypeT &resolver = *intercommunicandData->second.first;
              for (unsigned int sharedFieldIndex = 0; sharedFieldIndex <= sharedObject.Values().size();
                  sharedFieldIndex++)
              {
                Send(resolver.Fields()[sharedFieldIndex].first,
                     resolver.Fields()[sharedFieldIndex].second,
                     label,
                     *sharedObject.Values()[sharedFieldIndex]);
              }
            }
          }
          void Receive(const std::string & fieldLabel,
                       RuntimeType type,
                       const std::string objectLabel,
                       hemelb::multiscale::BaseSharedValue & value)
          {
            if (type == RuntimeTypeTraits::GetType<double>())
            {
              static_cast<hemelb::multiscale::SharedValue<double> &>(value) = doubleContents[objectLabel
                  + "_" + fieldLabel];

            }

          }
          void Send(const std::string & fieldLabel,
                    RuntimeType type,
                    const std::string objectLabel,
                    hemelb::multiscale::BaseSharedValue & value)
          {
            if (type == RuntimeTypeTraits::GetType<double>())
            {
              doubleContents[objectLabel + "_" + fieldLabel] =
                  static_cast<hemelb::multiscale::SharedValue<double> &>(value);
            }

          }
          std::map<std::string, double> &doubleContents;
          double currentTime;
      };
    }
  }
}

#endif // HEMELB_UNITTEST_MULTISCALE_MOCKINTERCOMMUNICATOR_H
