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
           * Read the TinyXML structure and set up the iolet, or write it to a TinyXML structure.
           *
           * @param iParent Parent XML element
           * @param iIsLoading Read or write?
           * @param simConfig The config object being read
           */
          void DoIO(TiXmlElement *iParent, bool iIsLoading, configuration::SimConfig* iSimConfig);

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
          double GetPressureGradientAmplitude() const;

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
          double GetPeriod() const;

          /**
           * Set the period of the zero average pressure gradient sine wave imposed
           *
           * @param per sinusoidal pressure profile period
           */
          void SetPeriod(const double& per);

          /**
           * Get the Womersley number characterising the pulsatile flow simulation
           *
           * @return Womersley number
           */
          double GetWomersleyNumber() const;

          /**
           * Set the Womersley number characterising the pulsatile flow simulation
           *
           * @param womNumber Womersley number
           */
          void SetWomersleyNumber(const double& womNumber);

        private:

          /**
           * Evaluates the Bessel function of the first kind order 0 at z. The solution is approximated
           * with numTerms terms of the series
           *
           *    J0(z) ~=  1 - 0.25*z*z/1!^2 + (0.25*z*z)^2/2!^2 - (0.25*z*z)^3/3!^2 + ...
           *
           * See formula 9.1.12 of Abramowitz and Stegun, Handbook of Mathematical Functions.
           *
           * At the time of writing the Boost implementation of the Bessel functions did not support
           * complex arguments.
           *
           * @param z point in the complex plane where we want to evaluate J0
           * @param numTerms number of terms in the series expansion used to approximate the function
           * @return J0(z)
           */
          std::complex<double> BesselJ0ComplexArgument(const std::complex<double>& z,
                                                       unsigned numTerms = 50) const;

          double pressureGradientAmplitude; ///< See class documentation
          double period; ///< See class documentation
          double womersleyNumber; ///< See class documentation
      };
    }
  }
}
#endif // HEMELB_LB_IOLETS_INOUTLETWOMERSLEYVELOCITY_H
