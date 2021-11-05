// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_TESTS_MULTISCALE_MOCKINTERCOMMUNICAND_H
#define HEMELB_TESTS_MULTISCALE_MOCKINTERCOMMUNICAND_H

#include "multiscale/Intercommunicand.h"
#include "multiscale/SharedValue.h"

namespace hemelb
{
  namespace tests
  {
    /***
     * A mock intercommunicand, representing an iolet which has density (a.k.a. pressure) and a velocity.
     */
    class MockIntercommunicand : public multiscale::Intercommunicand
    {
    public:
      MockIntercommunicand(double pressure=0,double velocity=0);
      virtual ~MockIntercommunicand();
      double GetPressure();
      double GetVelocity();
      void SetPressure(double aPressure);
      void SetVelocity(double aVelocity);

    private:
      multiscale::SharedValue<double> pressure;
      multiscale::SharedValue<double> velocity;

    };
  }
}

#endif //HEMELB_TESTS_MULTISCALE_MOCKINTERCOMMUNICAND_H
