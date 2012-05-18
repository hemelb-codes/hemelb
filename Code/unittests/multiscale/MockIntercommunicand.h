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
            pressure=aPressure;
          }
          void SetVelocity(double aVelocity){
            velocity=aVelocity;
          }
        private:
          SharedValue<double> pressure;
          SharedValue<double> velocity;

      };
    }
  }
}

#endif //HEMELB_UNITTESTS_MULTISCALE_MOCKINTERCOMMUNICAND_H
