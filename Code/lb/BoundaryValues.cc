#include "lb/BoundaryValues.h"
#include "topology/NetworkTopology.h"
#include "util/utilityFunctions.h"
#include <math.h>

namespace hemelb
{
  namespace lb
  {

    BoundaryValues::BoundaryValues(BoundaryComms* iInletComms,
                                   BoundaryComms* iOutletComms,
                                   geometry::LatticeData* iLatDat,
                                   SimConfig* iSimConfig,
                                   SimulationState* iSimState,
                                   util::UnitConverter* iUnits) :
      mState(iSimState), mSimConfig(iSimConfig), mUnits(iUnits)
    {
      nTotInlets = (int) iSimConfig->Inlets.size();
      nTotOutlets = (int) iSimConfig->Outlets.size();

      ReadParameters();

      InitialiseBoundaryDensities();

      iInletComms->Initialise(geometry::LatticeData::INLET_TYPE, iLatDat, inlet_density_cycle);

      iOutletComms->Initialise(geometry::LatticeData::OUTLET_TYPE, iLatDat, outlet_density_cycle);
    }

    BoundaryValues::~BoundaryValues()
    {

    }

    void BoundaryValues::InitialiseBoundaryDensities()
    {
      if (topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
      {
        double w = 2.0 * PI / (double) mState->GetTimeStepsPerCycle();

        for (unsigned long time_step = 0; time_step < mState->GetTimeStepsPerCycle(); time_step++)
        {
          for (int i = 0; i < nTotInlets; i++)
          {
            inlet_density_cycle[time_step * nTotInlets + i] = inlet_density_avg[i]
                + inlet_density_amp[i] * cos(w * (double) time_step + inlet_density_phs[i]);
          }
          for (int i = 0; i < nTotOutlets; i++)
          {
            outlet_density_cycle[time_step * nTotOutlets + i] = outlet_density_avg[i]
                + outlet_density_amp[i] * cos(w * (double) time_step + outlet_density_phs[i]);
          }
        }
      }
    }

    void BoundaryValues::ReadParameters()
    {
      allocateInlets();

      for (int n = 0; n < nTotInlets; n++)
      {
        hemelb::SimConfig::InOutLet *lInlet = &mSimConfig->Inlets[n];

        inlet_density_avg[n] = mUnits->ConvertPressureToLatticeUnits(lInlet->PMean) / Cs2;
        inlet_density_amp[n] = mUnits->ConvertPressureGradToLatticeUnits(lInlet->PAmp) / Cs2;
        inlet_density_phs[n] = lInlet->PPhase * DEG_TO_RAD;
      }

      allocateOutlets();

      for (int n = 0; n < nTotOutlets; n++)
      {
        hemelb::SimConfig::InOutLet *lOutlet = &mSimConfig->Outlets[n];
        outlet_density_avg[n] = mUnits->ConvertPressureToLatticeUnits(lOutlet->PMean) / Cs2;
        outlet_density_amp[n] = mUnits->ConvertPressureGradToLatticeUnits(lOutlet->PAmp) / Cs2;
        outlet_density_phs[n] = lOutlet->PPhase * DEG_TO_RAD;
      }

    }

    void BoundaryValues::allocateInlets()
    {
      if (topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
      {
        inlet_density_cycle = new distribn_t[util::NumericalFunctions::max<int>(1, nTotInlets)
            * mState->GetTimeStepsPerCycle()];
      }

      inlet_density_avg = new distribn_t[nTotInlets];
      inlet_density_amp = new distribn_t[nTotInlets];
      inlet_density_phs = new distribn_t[nTotInlets];
    }

    void BoundaryValues::allocateOutlets()
    {
      if (topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
      {
        outlet_density_cycle = new distribn_t[util::NumericalFunctions::max<int>(1, nTotOutlets)
            * mState->GetTimeStepsPerCycle()];
      }

      outlet_density_avg = new distribn_t[nTotOutlets];
      outlet_density_amp = new distribn_t[nTotOutlets];
      outlet_density_phs = new distribn_t[nTotOutlets];
    }

    distribn_t BoundaryValues::GetInitialDensity()
    {
      distribn_t density;

      density = 0.0;

      for (int i = 0; i < nTotOutlets; i++)
      {
        density += outlet_density_avg[i] - outlet_density_amp[i];
      }

      density /= nTotOutlets;

      return density;
    }

    distribn_t BoundaryValues::GetInletDensityMin(int iBoundaryId)
    {
      return inlet_density_avg[iBoundaryId] - inlet_density_amp[iBoundaryId];
    }

    distribn_t BoundaryValues::GetInletDensityMax(int iBoundaryId)
    {
      return inlet_density_avg[iBoundaryId] + inlet_density_amp[iBoundaryId];
    }

    distribn_t BoundaryValues::GetOutletDensityMin(int iBoundaryId)
    {
      return outlet_density_avg[iBoundaryId] - outlet_density_amp[iBoundaryId];
    }

    distribn_t BoundaryValues::GetOutletDensityMax(int iBoundaryId)
    {
      return outlet_density_avg[iBoundaryId] + outlet_density_amp[iBoundaryId];
    }

    void BoundaryValues::Reset()
    {
      int i;

      for (i = 0; i < nTotInlets; i++)
      {
        inlet_density_avg[i] = mUnits->ConvertPressureToPhysicalUnits(inlet_density_avg[i] * Cs2);
        inlet_density_amp[i] = mUnits->ConvertPressureGradToPhysicalUnits(inlet_density_amp[i]
            * Cs2);
      }
      for (i = 0; i < nTotOutlets; i++)
      {
        outlet_density_avg[i] = mUnits->ConvertPressureToPhysicalUnits(outlet_density_avg[i] * Cs2);
        outlet_density_amp[i] = mUnits->ConvertPressureGradToPhysicalUnits(outlet_density_amp[i]
            * Cs2);
      }

      mState->DoubleTimeResolution();

      for (i = 0; i < nTotInlets; i++)
      {
        inlet_density_avg[i] = mUnits->ConvertPressureToLatticeUnits(inlet_density_avg[i]) / Cs2;
        inlet_density_amp[i] = mUnits->ConvertPressureGradToLatticeUnits(inlet_density_amp[i])
            / Cs2;
      }
      for (i = 0; i < nTotOutlets; i++)
      {
        outlet_density_avg[i] = mUnits->ConvertPressureToLatticeUnits(outlet_density_avg[i]) / Cs2;
        outlet_density_amp[i] = mUnits->ConvertPressureGradToLatticeUnits(outlet_density_amp[i])
            / Cs2;
      }

      if (topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
      {
        delete[] inlet_density_cycle;
        delete[] outlet_density_cycle;

        inlet_density_cycle = new distribn_t[hemelb::util::NumericalFunctions::max<int>(1,
                                                                                        nTotInlets)
            * mState->GetTimeStepsPerCycle()];
        outlet_density_cycle
            = new distribn_t[hemelb::util::NumericalFunctions::max<int>(1, nTotOutlets)
                * mState->GetTimeStepsPerCycle()];
      }

      InitialiseBoundaryDensities();
    }

    // Calculate the BCs for each boundary site type and the
    // non-equilibrium distribution functions.
    void BoundaryValues::CalculateBC(distribn_t f[],
                                     hemelb::geometry::LatticeData::SiteType iSiteType,
                                     unsigned int iBoundaryId,
                                     distribn_t *density,
                                     distribn_t *vx,
                                     distribn_t *vy,
                                     distribn_t *vz,
                                     distribn_t f_neq[],
                                     BoundaryComms* iInletComms,
                                     BoundaryComms* iOutletComms)
    {
      distribn_t dummy_density;

      for (unsigned int l = 0; l < D3Q15::NUMVECTORS; l++)
      {
        f_neq[l] = f[l];
      }

      if (iSiteType == hemelb::geometry::LatticeData::FLUID_TYPE)
      {
        D3Q15::CalculateDensityAndVelocity(f, *density, *vx, *vy, *vz);
      }
      else
      {
        if (iSiteType == hemelb::geometry::LatticeData::INLET_TYPE)
        {
          *density = iInletComms->GetBoundaryDensity(iBoundaryId);
        }
        else
        {
          *density = iOutletComms->GetBoundaryDensity(iBoundaryId);
        }

        D3Q15::CalculateDensityAndVelocity(f, dummy_density, *vx, *vy, *vz);
        D3Q15::CalculateFeq(*density, *vx, *vy, *vz, f);

      }
      for (unsigned int l = 0; l < D3Q15::NUMVECTORS; l++)
      {
        f_neq[l] -= f[l];
      }

    }

  }
}
