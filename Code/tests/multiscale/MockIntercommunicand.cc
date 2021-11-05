// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "tests/multiscale/MockIntercommunicand.h"

namespace hemelb
{
  namespace tests
  {
    /***
     * A mock intercommunicand, representing an iolet which has density (a.k.a. pressure) and a velocity.
     */
    
    MockIntercommunicand::MockIntercommunicand(double pressure, double velocity) :
      Intercommunicand(), pressure(this,pressure), velocity(this,velocity)
    {
    }
    MockIntercommunicand::~MockIntercommunicand()
    {
    }
    double MockIntercommunicand::GetPressure(){
      return pressure;
    }
    double MockIntercommunicand::GetVelocity(){
      return velocity;
    }
    void MockIntercommunicand::SetPressure(double aPressure){
      pressure.SetPayload(aPressure);
    }
    void MockIntercommunicand::SetVelocity(double aVelocity){
      velocity.SetPayload(aVelocity);
    }
  }
}

