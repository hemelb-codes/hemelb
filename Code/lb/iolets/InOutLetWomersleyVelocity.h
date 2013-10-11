//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_LB_IOLETS_INOUTLETWOMERSLEYVELOCITY_H
#define HEMELB_LB_IOLETS_INOUTLETWOMERSLEYVELOCITY_H
#include "lb/iolets/InOutLetVelocity.h"
#include <complex>

namespace hemelb
{
  namespace lb
  {
    namespace iolets
    {
      /**
       * This class implements the InOutLetVelocity interface. It imposes a Womersley velocity profile
       * at a given inlet/outlet of a simulation of laminar flow in a cylinder of radius R and zero
       * average sinusoidal pressure gradient given by
       *
       *   pressureGradientAmplitude * sin(2*pi*t/period)
       *
       * where t is the current time. See physics validation paper for the analytical expression of
       * velocity implemented by GetVelocity(x, t) which is also a function of pressureGradientAmplitude,
       * radius, period, and womersleyNumber.
       *
       * If combined with a pressure iolet at the other end of the cylinder, it must be set to
       * zero pressure
       */
      class InOutLetWomersleyVelocity : public InOutLetVelocity
      {
        public:

          /**
           * Returns a copy of the current iolet. The caller is responsible for freeing that memory.
           *
           * @return copy of the current iolet
           */
          InOutLet* Clone() const;

          /**
           * Get Womersley velocity for a given time and position.
           *
           * @param x lattice site position
           * @param t time
           * @return velocity
           */
          LatticeVelocity GetVelocity(const LatticePosition& x, const LatticeTime t) const;

          /**
           * Get the amplitude of the zero average pressure gradient sine wave imposed.
           *
           * @return pressure gradient amplitude
           */
          double& GetPressureGradientAmplitude();

          /**
           * Set the amplitude of the zero average pressure gradient sine wave imposed.
           *
           * @param pressGradAmp pressure gradient amplitude
           */
          void SetPressureGradientAmplitude(const double& pressGradAmp);

          /**
           * Get the period of the zero average pressure gradient sine wave imposed
           *
           * @return sinusoidal pressure profile period
           */
          LatticeTime& GetPeriod();

          /**
           * Set the period of the zero average pressure gradient sine wave imposed
           *
           * @param per sinusoidal pressure profile period
           */
          void SetPeriod(const LatticeTime& per);

          /**
           * Get the Womersley number characterising the pulsatile flow simulation
           *
           * @return Womersley number
           */
          Dimensionless& GetWomersleyNumber();

          /**
           * Set the Womersley number characterising the pulsatile flow simulation
           *
           * @param womNumber Womersley number
           */
          void SetWomersleyNumber(const Dimensionless& womNumber);

        private:
          typedef std::complex<double> Complex;
          static const Complex i;
          static const Complex iPowThreeHalves;
          double pressureGradientAmplitude; ///< See class documentation
          LatticeTime period; ///< See class documentation
          double womersleyNumber; ///< See class documentation
      };
    }
  }
}
#endif // HEMELB_LB_IOLETS_INOUTLETWOMERSLEYVELOCITY_H
