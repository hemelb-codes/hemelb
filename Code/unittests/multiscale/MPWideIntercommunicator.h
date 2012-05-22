#ifndef HEMELB_UNITTESTS_MULTISCALE_MOCKINTERCOMMUNICATOR_H
#define HEMELB_UNITTESTS_MULTISCALE_MOCKINTERCOMMUNICATOR_H
#include "multiscale/Intercommunicator.h"
#include "mpiInclude.h"
#include "tinyxml.h"

/**
MPWide Intercommunicator class.
This class provides an abstraction for the MPWide Communication Library.
http://castle.strw.leidenuniv.nl/software/MPWide.html
**/

#include <algorithm>
#include <functional>
#include <vector>
#include <sstream>
#include "MPWide.h"

/*
TODO: Define BoundaryRegion - a unified structure for boundary regions in HemeLB.


*/


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
      class MPWideIntercommunicator : public hemelb::multiscale::Intercommunicator<MPIRuntimeType>
      {
        public:
          MPWideIntercommunicator(std::map<std::string, double> & buffer,std::map<std::string,bool> &orchestration) :
              doubleContents(buffer), currentTime(0), orchestration(orchestration)
          {
              //TODO: Decide where we are going to fetch this information (separate config file?)
              string *url
              int* server_side_ports
              int num_channels

              MPW_Init(url, server_side_ports, num_channels);
          }

          void ShareInitialConditions()            
          {                                        
            doubleContents["shared_time"] = 0.0;   
            SendToMultiscale();                    
          }                                        
                                         
          bool DoMultiscale(double new_time)
          {      
            if (ShouldAdvance()) {                  
              AdvanceTime(new_time);               
              SendToMultiscale();                  
            }                                      
            bool should_advance=ShouldAdvance();   
            if (should_advance) {                                      
              GetFromMultiscale();                 
            }                                      
            return should_advance;                 
          }                                        
        private:
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

              for (unsigned int sharedFieldIndex = 0; sharedFieldIndex < sharedObject.Values().size();
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
              for (unsigned int sharedFieldIndex = 0; sharedFieldIndex < sharedObject.Values().size();
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

            std::string label(objectLabel+ "_" + fieldLabel);
            if (orchestration[label]) return;
            if (type == RuntimeTypeTraits::GetType<double>())
            {  
              static_cast<hemelb::multiscale::SharedValue<double> &>(value) = doubleContents[label];
            }

            /* For latency reasons we want to receive full boundary information in one data type, 
               using two calls at most. This is why we introduce a BoundaryRegion datatype.
               Later versions may be able to merge this into one call using MPW_DSendRecv(). */
            if (type == RuntimeTypeTraits::GetType<BoundaryRegion>())
            {
              long long int br_size = 0;
              MPW_Recv((char *) &br_size, sizeof(long long int), channels, num_channels);

              char *br_serialized = (char *) malloc(br_size);

              MPW_Recv(br_serialized, br_size, channels, num_channels);

              DeSerializeBoundaryRegion(br_serialized, value);
            }

          }
          void Send(const std::string & fieldLabel,
                    RuntimeType type,
                    const std::string objectLabel,
                    hemelb::multiscale::BaseSharedValue & value)
          {
            std::string label(objectLabel+ "_" + fieldLabel);
            if (!orchestration[label]) return;
            if (type == RuntimeTypeTraits::GetType<double>())
            {
              doubleContents[label] =
                  static_cast<hemelb::multiscale::SharedValue<double> &>(value);
            }

            if (type == RuntimeTypeTraits::GetType<BoundaryRegion>())
            {
                long long int br_size = 0;
                MPW_Recv((char *) &br_size, sizeof(long long int), channels, num_channels);

                char *br_serialized = (char *) malloc(br_size);

                MPW_Recv(br_serialized, br_size, channels, num_channels);
  
                DeSerializeBoundaryRegion(br_serialized, value);
            }


          }
          std::map<std::string, double> &doubleContents;
          double currentTime;
          std::map<std::string,bool> & orchestration;
      };
    }
  }
}

#endif // HEMELB_UNITTEST_MULTISCALE_MOCKINTERCOMMUNICATOR_H
