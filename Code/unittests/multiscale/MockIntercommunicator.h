// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

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
       * Orchestration is true if it is intent out, and false if intent in.
       */
      class MockIntercommunicator : public hemelb::multiscale::Intercommunicator<MPIRuntimeType>
      {
        public:
          MockIntercommunicator(std::map<std::string, double> & buffer,std::map<std::string,bool> &orchestration) :
              doubleContents(buffer), currentTime(0), orchestration(orchestration)
          {

          }

          void ShareInitialConditions()
          {
            doubleContents["shared_time"] = 0.0;
            SendToMultiscale();
          }

          bool DoMultiscale(double new_time){
            if (ShouldAdvance()){
              AdvanceTime(new_time);
              SendToMultiscale();
            }
            bool should_advance=ShouldAdvance();
            if (should_advance)
            {
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

              for (unsigned int sharedFieldIndex = 0; sharedFieldIndex < sharedObject.SharedValues().size();
                  sharedFieldIndex++)
              {

                Receive(resolver.Fields()[sharedFieldIndex].first,
                        resolver.Fields()[sharedFieldIndex].second,
                        label,
                        *sharedObject.SharedValues()[sharedFieldIndex]);
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
              for (unsigned int sharedFieldIndex = 0; sharedFieldIndex < sharedObject.SharedValues().size();
                  sharedFieldIndex++)
              {
                Send(resolver.Fields()[sharedFieldIndex].first,
                     resolver.Fields()[sharedFieldIndex].second,
                     label,
                     *sharedObject.SharedValues()[sharedFieldIndex]);
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
              (static_cast<hemelb::multiscale::SharedValue<double> &>(value)).SetPayload(doubleContents[label]);
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

          }
          std::map<std::string, double> &doubleContents;
          double currentTime;
          std::map<std::string,bool> & orchestration;
      };
    }
  }
}

#endif // HEMELB_UNITTEST_MULTISCALE_MOCKINTERCOMMUNICATOR_H
