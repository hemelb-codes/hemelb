#ifndef HEMELB_LB_LBMCONFIG_H
#define HEMELB_LB_LBMCONFIG_H

namespace hemelb
{
  namespace lb
  {
    enum StressTypes
    {
      VonMises = 0,
      ShearStress = 1,
      IgnoreStress = 2
    };

    struct LbmParameters
    {
        double Omega;
        double Tau;
        double StressParameter;
        double Beta;  // Viscous dissipation in ELBM
        StressTypes StressType;
    };
  }
}

#endif //HEMELB_LB_LBMCONFIG_H
