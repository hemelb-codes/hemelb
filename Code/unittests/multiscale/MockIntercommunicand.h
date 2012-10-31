// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_UNITTESTS_MULTISCALE_MOCKINTERCOMMUNICAND_H
#define HEMELB_UNITTESTS_MULTISCALE_MOCKINTERCOMMUNICAND_H

#include "multiscale/Intercommunicand.h"
#include "multiscale/SharedValue.h"
#include <iostream>
#include <vector>

namespace hemelb
{
  namespace unittests
  {
    namespace multiscale
    {
      using namespace hemelb::multiscale;
      /***
       * A mock intercommunicand, representing an iolet which has density (a.k.a. pressure) and a velocity.
       */
      class MockIntercommunicand : public Intercommunicand
      {
        public:
          MockIntercommunicand(double pressure=0,double velocity=0) :
          Intercommunicand(), pressure(this,pressure), velocity(this,velocity)
          {
          }
          virtual ~MockIntercommunicand()
          {
          }
          ;
          double GetPressure(){
            return pressure;
          }
          double GetVelocity(){
            return velocity;
          }
          void SetPressure(double aPressure){
            pressure.SetPayload(aPressure);
          }
          void SetVelocity(double aVelocity){
            velocity.SetPayload(aVelocity);
          }
        private:
          SharedValue<double> pressure;
          SharedValue<double> velocity;

      };
    }
  }
}

#endif //HEMELB_UNITTESTS_MULTISCALE_MOCKINTERCOMMUNICAND_H
