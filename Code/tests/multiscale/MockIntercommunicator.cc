// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "tests/multiscale/MockIntercommunicator.h"

#include "multiscale/SharedValue.h"

namespace hemelb
{
  namespace tests
  {

    MockIntercommunicator::MockIntercommunicator(std::map<std::string, double> & buffer, std::map<std::string, bool> &orchestration) :
      doubleContents(buffer), currentTime(0), orchestration(orchestration)
    {

    }

    void MockIntercommunicator::ShareInitialConditions()
    {
      doubleContents["shared_time"] = 0.0;
      SendToMultiscale();
    }

    bool MockIntercommunicator::DoMultiscale(double new_time)
    {
      if (ShouldAdvance())
	{
	  AdvanceTime(new_time);
	  SendToMultiscale();
	}
      bool should_advance = ShouldAdvance();
      if (should_advance)
	{
	  GetFromMultiscale();
	}
      return should_advance;
    }
    
    void MockIntercommunicator::AdvanceTime(double new_time)
    {
      currentTime = new_time;
      doubleContents["shared_time"] = new_time;
    }
    bool MockIntercommunicator::ShouldAdvance()
    {
      return doubleContents["shared_time"] >= currentTime;
    }

    bool MockIntercommunicator::GetFromMultiscale()
    {
      for (ContentsType::iterator intercommunicandData = registeredObjects.begin();
	   intercommunicandData != registeredObjects.end(); intercommunicandData++)
	{
	  multiscale::Intercommunicand &sharedObject = *intercommunicandData->first;
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
    void MockIntercommunicator::SendToMultiscale()
    {
      for (ContentsType::iterator intercommunicandData = registeredObjects.begin();
	   intercommunicandData != registeredObjects.end(); intercommunicandData++)
	{
	  multiscale::Intercommunicand &sharedObject = *intercommunicandData->first;
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
    void MockIntercommunicator::Receive(const std::string & fieldLabel,
					RuntimeType type,
					const std::string objectLabel,
					multiscale::BaseSharedValue & value)
    {

      std::string label(objectLabel + "_" + fieldLabel);
      if (orchestration[label])
	return;
      if (type == RuntimeTypeTraits::GetType<double>())
	{
	  (static_cast<multiscale::SharedValue<double> &>(value)).SetPayload(doubleContents[label]);
	}

    }
    void MockIntercommunicator::Send(const std::string & fieldLabel,
				     RuntimeType type,
				     const std::string objectLabel,
				     multiscale::BaseSharedValue & value)
    {
      std::string label(objectLabel + "_" + fieldLabel);
      if (!orchestration[label])
	return;
      if (type == RuntimeTypeTraits::GetType<double>())
	{
	  doubleContents[label] = static_cast<multiscale::SharedValue<double> &>(value);
	}

    }

  }
}

