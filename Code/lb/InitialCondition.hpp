// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_INITIALCONDITION_HPP
#define HEMELB_LB_INITIALCONDITION_HPP


namespace hemelb {
  namespace lb {
    template<class LatticeType>
    struct FSetter {
      using result_type = void;
      template <typename T>
      void operator()(T t) const {
	t.template SetFs<LatticeType>(latDat, ioComms);
      }
      geometry::LatticeData* latDat;
      const net::IOCommunicator& ioComms;
    };

    
    template<class LatticeType>
    void InitialCondition::SetFs(geometry::LatticeData* latDat, const net::IOCommunicator& ioComms) const {
      boost::apply_visitor(FSetter<LatticeType>{latDat, ioComms}, *this);
    }

    template<class LatticeType>
    void EquilibriumInitialCondition::SetFs(geometry::LatticeData* latDat, const net::IOCommunicator& ioComms) const {
      distribn_t f_eq[LatticeType::NUMVECTORS];
      LatticeType::CalculateFeq(density, mom_x, mom_y, mom_z, f_eq);
      
      for (site_t i = 0; i < latDat->GetLocalFluidSiteCount(); i++) {
	distribn_t* f_old_p = this->GetFOld(latDat, i * LatticeType::NUMVECTORS);
	distribn_t* f_new_p = this->GetFNew(latDat, i * LatticeType::NUMVECTORS);
	
	for (unsigned int l = 0; l < LatticeType::NUMVECTORS; l++) {
	  f_new_p[l] = f_old_p[l] = f_eq[l];
	}
      }
    }


    template<class LatticeType>
    void CheckpointInitialCondition::SetFs(geometry::LatticeData* latDat, const net::IOCommunicator& ioComms) const {
      auto distributionInputPtr = std::make_unique<extraction::LocalDistributionInput>(cpFile, ioComms);
      distributionInputPtr->LoadDistribution(latDat, initial_time);
    }

  }
}

#endif
