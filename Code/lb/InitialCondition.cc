// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "lb/InitialCondition.h"

namespace hemelb {
  namespace lb {

    // InitialConditionBase
    
    InitialConditionBase::InitialConditionBase() {
    }
    InitialConditionBase::InitialConditionBase(boost::optional<LatticeTimeStep> t) : initial_time(t) {
    }
    
    void InitialConditionBase::SetTime(SimulationState* sim) const {
      if (initial_time)
	sim->timeStep = *initial_time;
    }


    // Equilibrium
    
    EquilibriumInitialCondition::EquilibriumInitialCondition() :
      InitialConditionBase(),
      density(1.0), mom_x(0.0), mom_y(0.0), mom_z(0.0) {
    }
    
    EquilibriumInitialCondition::EquilibriumInitialCondition(
      boost::optional<LatticeTimeStep> t0,
      distribn_t rho,
      distribn_t mx, distribn_t my, distribn_t mz) :
      InitialConditionBase(t0),
      density(rho),
      mom_x(mx), mom_y(my), mom_z(mz) {
    }
    
    CheckpointInitialCondition::CheckpointInitialCondition(boost::optional<LatticeTimeStep> t0, const std::string& cp)
      : InitialConditionBase(t0), cpFile(cp) {
    }

    // InitialCondition - sum type container

    // Visitor for setting time
    struct TSetter {
      using result_type = void;
      template <typename T>
      void operator()(T t) const {
	t.SetTime(ss);
      }
      SimulationState* ss;
    };
    void InitialCondition::SetTime(SimulationState* sim) const {
      const ICVar* self = this;
      boost::apply_visitor(TSetter{sim}, *self);
    }

    // See InitialCondtions.hpp for setting Fs (distributions)

    // Visitor for factory function
    struct ICMaker {
      using result_type = InitialCondition;
      
      template <typename T>
      InitialCondition operator()(T) const {
	throw Exception() << "Trying to make an InitialCondition from unknown type of config";
      }

      InitialCondition operator()(const configuration::EquilibriumIC& cfg) const {
	auto rho = cfg.unitConverter->ConvertPressureToLatticeUnits(cfg.p_mmHg) / Cs2;
	return EquilibriumInitialCondition{cfg.t0, rho};
      }
      InitialCondition operator()(const configuration::CheckpointIC& cfg) const {
	return CheckpointInitialCondition{cfg.t0, cfg.cpFile};
      }
    };
    
    // Factory function just delegates to visitor
    InitialCondition InitialCondition::FromConfig(const configuration::ICConfig& conf) {
      return boost::apply_visitor(ICMaker{}, conf);
    }
  }
}
