#include "lb/BoundaryValues.h"
#include "topology/NetworkTopology.h"
#include "util/utilityFunctions.h"
#include "util/fileutils.h"
#include <fstream>
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

      inlet_density_cycle = std::vector<distribn_t>(0);
      outlet_density_cycle = std::vector<distribn_t>(0);

      ReadParameters();

      InitialiseBoundaryDensities();

      iInletComms->Initialise(geometry::LatticeData::INLET_TYPE, iLatDat, &inlet_density_cycle);

      iOutletComms->Initialise(geometry::LatticeData::OUTLET_TYPE, iLatDat, &outlet_density_cycle);
    }

    BoundaryValues::~BoundaryValues()
    {
      delete[] inlet_density_avg;
      delete[] inlet_density_amp;
      delete[] inlet_density_phs;
      delete[] inlet_density_min;
      delete[] inlet_density_max;
      delete[] inlet_file;

      delete[] outlet_density_avg;
      delete[] outlet_density_amp;
      delete[] outlet_density_phs;
      delete[] outlet_density_min;
      delete[] outlet_density_max;
      delete[] outlet_file;
    }

    void BoundaryValues::InitialiseBoundaryDensities()
    {
      if (topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
      {

        for (int i = 0; i < nTotInlets; i++)
        {
          if (inlet_file[i] == "")
          {
            InitialiseCosCycle(i,
                               nTotInlets,
                               inlet_density_avg,
                               inlet_density_amp,
                               inlet_density_phs,
                               inlet_density_cycle);
          }
          else
          {
            InitialiseFromFile(i, inlet_file[i], inlet_density_cycle);
          }
        }

        for (int i = 0; i < nTotOutlets; i++)
        {
          if (outlet_file[i] == "")
          {
            InitialiseCosCycle(i,
                               nTotOutlets,
                               outlet_density_avg,
                               outlet_density_amp,
                               outlet_density_phs,
                               outlet_density_cycle);
          }
          else
          {
            InitialiseFromFile(i, outlet_file[i], outlet_density_cycle);
          }
        }

        FindIOletDensityExtrema();
      }

      MPI_Bcast(inlet_density_min,
                nTotInlets,
                hemelb::MpiDataType(inlet_density_min[0]),
                0,
                MPI_COMM_WORLD);

      MPI_Bcast(inlet_density_max,
                nTotInlets,
                hemelb::MpiDataType(inlet_density_max[0]),
                0,
                MPI_COMM_WORLD);

      MPI_Bcast(outlet_density_min,
                nTotOutlets,
                hemelb::MpiDataType(outlet_density_min[0]),
                0,
                MPI_COMM_WORLD);

      MPI_Bcast(outlet_density_max,
                nTotOutlets,
                hemelb::MpiDataType(outlet_density_max[0]),
                0,
                MPI_COMM_WORLD);
    }

    void BoundaryValues::InitialiseCosCycle(int i,
                                            int IOlets,
                                            distribn_t* density_avg,
                                            distribn_t* density_amp,
                                            distribn_t* density_phs,
                                            std::vector<distribn_t> &density_cycle)
    {
      double w = 2.0 * PI / (double) mState->GetTimeStepsPerCycle();

      for (unsigned long time_step = 0; time_step < mState->GetTimeStepsPerCycle(); time_step++)
      {
        density_cycle[time_step * IOlets + i] = density_avg[i] + density_amp[i] * cos(w
            * (double) time_step + density_phs[i]);
      }
    }

    void BoundaryValues::InitialiseFromFile(int i,
                                            std::string &filename,
                                            std::vector<distribn_t> &density_cycle)
    {
      // First read in values from file into vectors

      std::vector<double> time(0);
      std::vector<double> value(0);

      double timeTemp, valueTemp;

      util::check_file(filename.c_str());
      std::ifstream datafile(filename.c_str());

      while (datafile.good())
      {
        datafile >> timeTemp >> valueTemp;
        time.push_back(timeTemp);
        value.push_back(valueTemp);
      }

      datafile.close();

      // Now convert these vectors into arrays using linear interpolation

      for (unsigned long time_step = 0; time_step < mState->GetTimeStepsPerCycle(); time_step++)
      {
        double point = time[0] + (double) time_step / (double) mState->GetTimeStepsPerCycle()
            * (time[time.size() - 1] - time[0]);

        density_cycle[time_step * nTotInlets + i]
            = mUnits->ConvertPressureToLatticeUnits(util::NumericalFunctions::LinearInterpolate(time,
                                                                                                value,
                                                                                                point))
                / Cs2;
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

        inlet_file[n] = lInlet->PFilePath;
      }

      allocateOutlets();

      for (int n = 0; n < nTotOutlets; n++)
      {
        hemelb::SimConfig::InOutLet *lOutlet = &mSimConfig->Outlets[n];

        outlet_density_avg[n] = mUnits->ConvertPressureToLatticeUnits(lOutlet->PMean) / Cs2;
        outlet_density_amp[n] = mUnits->ConvertPressureGradToLatticeUnits(lOutlet->PAmp) / Cs2;
        outlet_density_phs[n] = lOutlet->PPhase * DEG_TO_RAD;

        outlet_file[n] = lOutlet->PFilePath;
      }

    }

    void BoundaryValues::allocateInlets()
    {
      if (topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
      {
        inlet_density_cycle.resize(util::NumericalFunctions::max<int>(1, nTotInlets)
            * mState->GetTimeStepsPerCycle());
      }

      inlet_density_avg = new distribn_t[nTotInlets];
      inlet_density_amp = new distribn_t[nTotInlets];
      inlet_density_phs = new distribn_t[nTotInlets];
      inlet_density_min = new distribn_t[nTotInlets];
      inlet_density_max = new distribn_t[nTotInlets];

      inlet_file = new std::string[nTotInlets];
    }

    void BoundaryValues::allocateOutlets()
    {
      if (topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
      {
        outlet_density_cycle.resize(util::NumericalFunctions::max<int>(1, nTotOutlets)
            * mState->GetTimeStepsPerCycle());
      }

      outlet_density_avg = new distribn_t[nTotOutlets];
      outlet_density_amp = new distribn_t[nTotOutlets];
      outlet_density_phs = new distribn_t[nTotOutlets];
      outlet_density_min = new distribn_t[nTotOutlets];
      outlet_density_max = new distribn_t[nTotOutlets];

      outlet_file = new std::string[nTotOutlets];
    }

    void BoundaryValues::FindIOletDensityExtrema()
    {
      for (int i = 0; i < nTotInlets; i++)
      {
        if (inlet_file[i] == "")
        {
          inlet_density_min[i] = inlet_density_cycle[0];
          inlet_density_max[i] = inlet_density_cycle[0];

          for (unsigned int j = 0; j < inlet_density_cycle.size(); j++)
          {
            inlet_density_min[i] = util::NumericalFunctions::min(inlet_density_min[i],
                                                                 inlet_density_cycle[j]);
            inlet_density_max[i] = util::NumericalFunctions::max(inlet_density_max[i],
                                                                 inlet_density_cycle[j]);
          }
        }
        else
        {
          inlet_density_min[i] = inlet_density_avg[i] - inlet_density_amp[i];
          inlet_density_max[i] = inlet_density_avg[i] + inlet_density_amp[i];
        }
      }
      for (int i = 0; i < nTotOutlets; i++)
      {
        if (outlet_file[i] == "")
        {
          outlet_density_min[i] = outlet_density_cycle[0];
          outlet_density_max[i] = outlet_density_cycle[0];

          for (unsigned int j = 0; j < outlet_density_cycle.size(); j++)
          {
            outlet_density_min[i] = util::NumericalFunctions::min(outlet_density_min[i],
                                                                  outlet_density_cycle[j]);
            outlet_density_max[i] = util::NumericalFunctions::max(outlet_density_max[i],
                                                                  outlet_density_cycle[j]);
          }
        }
        else
        {
          outlet_density_min[i] = outlet_density_avg[i] - outlet_density_amp[i];
          outlet_density_max[i] = outlet_density_avg[i] + outlet_density_amp[i];
        }
      }
    }

    distribn_t BoundaryValues::GetInitialDensity()
    {
      distribn_t density;

      density = 0.0;

      for (int i = 0; i < nTotOutlets; i++)
      {
        density += GetOutletDensityMin(i);
      }

      density /= nTotOutlets;

      return density;
    }

    distribn_t BoundaryValues::GetInletDensityMin(int iBoundaryId)
    {
      return inlet_density_min[iBoundaryId];
    }

    distribn_t BoundaryValues::GetInletDensityMax(int iBoundaryId)
    {
      return inlet_density_max[iBoundaryId];
    }

    distribn_t BoundaryValues::GetOutletDensityMin(int iBoundaryId)
    {
      return outlet_density_min[iBoundaryId];
    }

    distribn_t BoundaryValues::GetOutletDensityMax(int iBoundaryId)
    {
      return outlet_density_max[iBoundaryId];
    }

    void BoundaryValues::Reset()
    {
      int i;

      for (i = 0; i < nTotInlets; i++)
      {
        if (inlet_file[i] != "")
        {
          inlet_density_avg[i] = mUnits->ConvertPressureToPhysicalUnits(inlet_density_avg[i] * Cs2);
          inlet_density_amp[i] = mUnits->ConvertPressureGradToPhysicalUnits(inlet_density_amp[i]
              * Cs2);
        }
      }
      for (i = 0; i < nTotOutlets; i++)
      {
        if (outlet_file[i] != "")
        {
          outlet_density_avg[i] = mUnits->ConvertPressureToPhysicalUnits(outlet_density_avg[i]
              * Cs2);
          outlet_density_amp[i] = mUnits->ConvertPressureGradToPhysicalUnits(outlet_density_amp[i]
              * Cs2);
        }
      }

      mState->DoubleTimeResolution();

      for (i = 0; i < nTotInlets; i++)
      {
        if (inlet_file[i] != "")
        {
          inlet_density_avg[i] = mUnits->ConvertPressureToLatticeUnits(inlet_density_avg[i]) / Cs2;
          inlet_density_amp[i] = mUnits->ConvertPressureGradToLatticeUnits(inlet_density_amp[i])
              / Cs2;
        }
      }
      for (i = 0; i < nTotOutlets; i++)
      {
        if (outlet_file[i] != "")
        {
          outlet_density_avg[i] = mUnits->ConvertPressureToLatticeUnits(outlet_density_avg[i])
              / Cs2;
          outlet_density_amp[i] = mUnits->ConvertPressureGradToLatticeUnits(outlet_density_amp[i])
              / Cs2;
        }
      }

      if (topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
      {
        inlet_density_cycle.resize(hemelb::util::NumericalFunctions::max<int>(1, nTotInlets)
            * mState->GetTimeStepsPerCycle());
        outlet_density_cycle.resize(hemelb::util::NumericalFunctions::max<int>(1, nTotOutlets)
            * mState->GetTimeStepsPerCycle());
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
