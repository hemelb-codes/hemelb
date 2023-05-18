// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_KERNELS_CARREAUYASUDARHEOLOGYMODEL_H
#define HEMELB_LB_KERNELS_CARREAUYASUDARHEOLOGYMODEL_H

#include "lb/kernels/AbstractRheologyModel.h"

namespace hemelb::lb
{
    struct InitParams;

    // Hold parameters defining a fit of the Carreau-Yasuda model.
    // Needs to be literal so it can be NTTP
    struct CYFit {
        double ETA_INF; // Pa.s
        double ETA_ZERO; // Pa.s
        double LAMBDA;// s
        double A; // dimensionless
        double N; // dimensionless
    };

    // Take the fit parameters as a NTTP to give compiler best chance to inline things.
    template<CYFit CYFIT>
    class CarreauYasudaRheologyModel : public AbstractRheologyModel<
            CarreauYasudaRheologyModel<CYFIT> >
    {
    public:
        CarreauYasudaRheologyModel(const InitParams&) {}
        /*
         *  Compute nu for a given shear rate according to the Carreau-Yasuda model:
         *
         *  eta = ETA_INF + (ETA_ZERO - ETA_INF) * (1 + (LAMBDA*iShearRate)^A)^((N-1)/A)
         *  nu = eta / density
         *
         *  @param iShearRate local shear rate value (s^{-1}).
         *  @param iDensity local density. TODO at the moment this value is not used
         *         in any subclass.
         *
         *  @return dynamic viscosity (Pa s).
         */
        double CalculateViscosityForShearRate(const double &iShearRate,
                                              const distribn_t &iDensity) const;
    };

    constexpr auto HumanCYFit = CYFit{0.0035, 0.16, 8.2, 0.64, 0.2128};
    constexpr auto MouseCYFit = CYFit{3.265e-3, 14.49e-3, 0.1829, 2.707, 0.4136};
    using CarreauYasudaRheologyModelHumanFit = CarreauYasudaRheologyModel<HumanCYFit>;
    using CarreauYasudaRheologyModelMouseFit = CarreauYasudaRheologyModel<MouseCYFit>;
}

#endif /* HEMELB_LB_KERNELS_CARREAUYASUDARHEOLOGYMODEL_H */
