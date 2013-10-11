#include "lb/iolets/InOutLetWomersleyVelocity.h"
#include "configuration/SimConfig.h"
#include "util/Bessel.h"

namespace hemelb
{
  namespace lb
  {
    namespace iolets
    {
      const InOutLetWomersleyVelocity::Complex InOutLetWomersleyVelocity::i = Complex(0, 1);
      const InOutLetWomersleyVelocity::Complex InOutLetWomersleyVelocity::iPowThreeHalves =
          pow(i, 1.5);

      InOutLet* InOutLetWomersleyVelocity::Clone() const
      {
        InOutLet* copy = new InOutLetWomersleyVelocity(*this);
        return copy;
      }

      LatticeVelocity InOutLetWomersleyVelocity::GetVelocity(const LatticePosition& x,
                                                             const LatticeTime t) const
      {
        LatticePosition posLat = units->ConvertPositionToLatticeUnits(position);
        LatticePosition displ = x - posLat;
        LatticeDistance z = displ.Dot(normal);
        Dimensionless r = sqrt(displ.GetMagnitudeSquared() - z * z);

        double omega = 2.0 * PI / period;
        LatticeDensity density = 1.0;

        Complex besselNumer = util::BesselJ0ComplexArgument(iPowThreeHalves * womersleyNumber * r
            / radius);
        Complex besselDenom = util::BesselJ0ComplexArgument(iPowThreeHalves * womersleyNumber);
        LatticeSpeed velocityMagnitude = std::real(pressureGradientAmplitude / (density * omega)
            * (1.0 - besselNumer / besselDenom) * exp(i * omega * double(t)));

        return normal * -velocityMagnitude;
      }

      double& InOutLetWomersleyVelocity::GetPressureGradientAmplitude()
      {
        return pressureGradientAmplitude;
      }

      void InOutLetWomersleyVelocity::SetPressureGradientAmplitude(const double& pressGradAmp)
      {
        pressureGradientAmplitude = pressGradAmp;
      }

      LatticeTime& InOutLetWomersleyVelocity::GetPeriod()
      {
        return period;
      }

      void InOutLetWomersleyVelocity::SetPeriod(const LatticeTime& per)
      {
        period = per;
      }

      Dimensionless& InOutLetWomersleyVelocity::GetWomersleyNumber()
      {
        return womersleyNumber;
      }

      void InOutLetWomersleyVelocity::SetWomersleyNumber(const Dimensionless& womNumber)
      {
        womersleyNumber = womNumber;
      }
    }
  }
}
