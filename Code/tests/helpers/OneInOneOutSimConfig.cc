// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "tests/helpers/OneInOneOutSimConfig.h"

namespace hemelb::tests::helpers
{
    OneInOneOutSimConfig::OneInOneOutSimConfig()
            : configuration::SimConfig(std::filesystem::path())
    {
        sim_info.time = {10000, 0, 60.0 / (70.0 * 1000.0)};
        sim_info.space = {0.01, util::Vector3D<PhysicalDistance>::Zero()};
        sim_info.fluid = {DEFAULT_FLUID_DENSITY_Kg_per_m3, DEFAULT_FLUID_VISCOSITY_Pas, 80.0};

        configuration::CosinePressureIoletConfig inlet;
        inlet.amp_mmHg = 1.0;
        inlet.mean_mmHg = 80.0;
        inlet.phase_rad = PI;
        inlet.period_s = 60.0 / 70.0;
        inlet.normal = util::Vector3D<Dimensionless>(-3, 4, -9);
        inlets.push_back(inlet);

        configuration::CosinePressureIoletConfig outlet;
        //auto outlet = util::make_clone_ptr<lb::iolets::InOutLetCosine>();
        outlet.amp_mmHg = 0.0;
        outlet.mean_mmHg = 80.0;
        outlet.phase_rad = 0.0;
        outlet.period_s = 60.0 / 70.0;
        outlet.normal = util::Vector3D<Dimensionless>(2, -1, 4);
        outlets.push_back(outlet);
    }

    void OneInOneOutSimConfig::CheckIoletMatchesCMake(const io::xml::Element& ioletEl,
                                                      const std::string& requiredBC) const
    {
    }

}
