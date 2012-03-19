#ifndef HEMELB_UNITTESTS_MULTISCALE_MOCKINTERCOMMUNICATOR_H
#define HEMELB_UNITTESTS_MULTISCALE_MOCKINTERCOMMUNICATOR_H
#include "multiscale/Intercommunicator.h"
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
      template<class T> struct ExampleTypeTraits
      {
          static unsigned int type;
      };

      template<class T> unsigned int ExampleTypeTraits<T>::type = 0;
      template<> unsigned int ExampleTypeTraits<double>::type = 1;
      template<> unsigned int ExampleTypeTraits<int>::type = 2;

      struct ExampleRuntimeTypeImplementation
      {
          typedef unsigned int RuntimeType;
          template<class T> static unsigned int GetType()
          {
            return ExampleTypeTraits<T>::type;
          }
      };

      class MockIntercommunicator : public hemelb::multiscale::Intercommunicator<ExampleRuntimeTypeImplementation>
      {
        public:
          MockIntercommunicator(std::map<std::string, double> & buffer) :
              double_contents(buffer), current_time(0)
          {

          }

          void AdvanceTime(double new_time)
          {
            current_time = new_time;
            double_contents["shared_time"] = new_time;
          }
          bool ShouldAdvance()
          {
            return double_contents["shared_time"] >= current_time;
          }

          bool GetFromMultiscale()
          {
            for (ContentsType::iterator intercommunicand_data = registered_objects.begin();
                intercommunicand_data != registered_objects.end(); intercommunicand_data++)
            {
              hemelb::multiscale::Intercommunicand &shared_object = *intercommunicand_data->first;
              std::string &label = intercommunicand_data->second.second;
              IntercommunicandTypeT &resolver = *intercommunicand_data->second.first;

              for (unsigned int shared_field_index = 0; shared_field_index <= shared_object.Values().size();
                  shared_field_index++)
              {
                Receive(resolver.Fields()[shared_field_index].first,
                        resolver.Fields()[shared_field_index].second,
                        label,
                        *shared_object.Values()[shared_field_index]);
              }
            }
            return true;
          }
          void SendToMultiscale()
          {
            for (ContentsType::iterator intercommunicand_data = registered_objects.begin();
                intercommunicand_data != registered_objects.end(); intercommunicand_data++)
            {
              hemelb::multiscale::Intercommunicand &shared_object = *intercommunicand_data->first;
              std::string &label = intercommunicand_data->second.second;
              IntercommunicandTypeT &resolver = *intercommunicand_data->second.first;
              for (unsigned int shared_field_index = 0; shared_field_index <= shared_object.Values().size();
                  shared_field_index++)
              {
                Send(resolver.Fields()[shared_field_index].first,
                     resolver.Fields()[shared_field_index].second,
                     label,
                     *shared_object.Values()[shared_field_index]);
              }
            }
          }
          void Receive(const std::string & field_label,
                       RuntimeType type,
                       const std::string object_label,
                       hemelb::multiscale::BaseSharedValue & value)
          {
            if (type == ExampleTypeTraits<double>::type)
            {
              static_cast<hemelb::multiscale::SharedValue<double> &>(value) = double_contents[object_label
                  + "_" + field_label];

            }

          }
          void Send(const std::string & field_label,
                    RuntimeType type,
                    const std::string object_label,
                    hemelb::multiscale::BaseSharedValue & value)
          {
            if (type == ExampleTypeTraits<double>::type)
            {
              double_contents[object_label + "_" + field_label] =
                  static_cast<hemelb::multiscale::SharedValue<double> &>(value);
            }

          }
          std::map<std::string, double> &double_contents;
          double current_time;
          double shared_time;
      };
    }
  }
}

#endif // HEMELB_UNITTEST_MULTISCALE_MOCKINTERCOMMUNICATOR_H
