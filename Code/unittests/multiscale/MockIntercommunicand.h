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
          MockIntercommunicand() :
          Intercommunicand(), pressure(this), velocity(this)
          {
          }
          virtual ~MockIntercommunicand()
          {
          }
          ;
          double GetDensity(){
            return pressure;
          }
          double GetVelocity(){
            return velocity;
          }
          void SetDensity(double aDensity){
            pressure=aDensity;
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
