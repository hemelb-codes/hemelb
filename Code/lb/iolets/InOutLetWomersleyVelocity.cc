#include "lb/iolets/InOutLetWomersleyVelocity.h"
#include "configuration/SimConfig.h"
#include <boost/math/special_functions/factorials.hpp>

namespace hemelb
{
  namespace lb
  {
    namespace iolets
    {
      void InOutLetWomersleyVelocity::DoIO(TiXmlElement *iParent, bool iIsLoading,
                                           configuration::SimConfig* iSimConfig)
      {
        iSimConfig->DoIOForWomersleyVelocityInOutlet(iParent, iIsLoading, this);
      }

      InOutLet* InOutLetWomersleyVelocity::Clone() const
      {
        InOutLet* copy = new InOutLetWomersleyVelocity(*this);
        return copy;
      }

      std::complex<double> InOutLetWomersleyVelocity::BesselJ0ComplexArgument(
          const std::complex<double>& z, unsigned numTerms) const
      {
        std::complex<double> sum(0.0, 0.0);

        for (unsigned i = 0; i < numTerms; i++)
        {
          double sign = pow(-1, i);
          sum += sign * pow(0.25 * z * z, i) / pow(boost::math::factorial<double>(i), 2);
        }
        // If this assertion trips, it is very likely that the power of two of the factorial overflowed
        assert( !std::isnan(real(sum)) && !std::isnan(imag(sum)));

        return sum;
      }

      LatticeVelocity InOutLetWomersleyVelocity::GetVelocity(const LatticePosition& x,
                                                             const LatticeTime t) const
      {
        LatticePosition posLat = units->ConvertPositionToLatticeUnits(position);
        LatticePosition displ = x - posLat;
        LatticeDistance z = displ.Dot(normal);
        Dimensionless r = sqrt(displ.GetMagnitudeSquared() - z * z);
        assert(r <= radius);

        double omega = 2.0 * PI / period;
        std::complex<double> i(0, 1);
        LatticeDensity density = 1.0;

        std::complex<double> besselNumer = BesselJ0ComplexArgument(pow(i, 1.5) * womersleyNumber * r
            / radius);
        std::complex<double> besselDenom = BesselJ0ComplexArgument(pow(i, 1.5) * womersleyNumber);
        LatticeSpeed velocityMagnitude = std::real(pressureGradientAmplitude / (density * omega)
            * (1.0 - besselNumer / besselDenom) * exp(i * omega * double(t)));

        return normal * -velocityMagnitude;
      }

      double InOutLetWomersleyVelocity::GetPressureGradientAmplitude() const
      {
        return pressureGradientAmplitude;
      }

      void InOutLetWomersleyVelocity::SetPressureGradientAmplitude(const double& pressGradAmp)
      {
        pressureGradientAmplitude = pressGradAmp;
      }

      double InOutLetWomersleyVelocity::GetPeriod() const
      {
        return period;
      }

      void InOutLetWomersleyVelocity::SetPeriod(const double& per)
      {
        period = per;
      }

      double InOutLetWomersleyVelocity::GetWomersleyNumber() const
      {
        return womersleyNumber;
      }

      void InOutLetWomersleyVelocity::SetWomersleyNumber(const double& womNumber)
      {
        womersleyNumber = womNumber;
      }
    }
  }
}
