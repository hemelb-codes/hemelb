#ifndef HEMELB_UNITTESTS_MULTISCALE_EXAMPLEINTERCOMMUNICAND_H
#define HEMELB_UNITTESTS_MULTISCALE_EXAMPLEINTERCOMMUNICAND_H

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
          void SetDensity(double adensity){
            pressure=adensity;
          }
          void SetVelocity(double avelocity){
            velocity=avelocity;
          }
        private:
          SharedValue<double> pressure;
          SharedValue<double> velocity;

      };
    }
  }
}

#endif //HEMELB_UNITTESTS_MULTISCALE_EXAMPLEINTERCOMMUNICAND_H
