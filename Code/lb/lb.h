#ifndef HEMELB_LB_H
#define HEMELB_LB_H

#include "net/net.h"
#include "topology/NetworkTopology.h"
#include "lb/collisions/Collisions.h"
#include "vis/ColPixel.h"
#include "SimConfig.h"

namespace hemelb
{
  namespace lb
  {
    class LBM
    {
      public:
        int total_fluid_sites;
        int inlets, outlets;
        int steering_session_id;
        int period;

        double ConvertPressureToLatticeUnits(double pressure) const;
        double ConvertPressureToPhysicalUnits(double pressure) const;
        double ConvertVelocityToLatticeUnits(double velocity) const;
        double ConvertStressToLatticeUnits(double stress) const;
        double ConvertStressToPhysicalUnits(double stress) const;
        double ConvertVelocityToPhysicalUnits(double velocity) const;

        LBM(hemelb::SimConfig *iSimulationConfig,
            const hemelb::topology::NetworkTopology * iNetTop,
            hemelb::lb::GlobalLatticeData &bGlobLatDat,
            int iSteeringSessionId,
            double* oFileReadTime);
        void Restart(hemelb::lb::LocalLatticeData &iLocalLatDat);
        ~LBM();

        int IsUnstable(hemelb::lb::LocalLatticeData &iLocalLatDat);

        hemelb::lb::Stability
        DoCycle(int perform_rt,
                net::Net *net,
                lb::LocalLatticeData *bLocallatDat,
                double &bLbTime,
                double &bMPISendTime,
                double &bMPIWaitTime);
        void CalculateFlowFieldValues();
        void RecalculateTauViscosityOmega();
        void UpdateBoundaryDensities(int cycle_id, int time_step);
        void
        UpdateInletVelocities(int time_step, lb::LocalLatticeData &iLocalLatDat, net::Net *net);

        void SetFTranslator(int* iFTranslator);

        void SetInitialConditions(hemelb::lb::LocalLatticeData &bLocalLatDat);

        void
        WriteConfig(hemelb::lb::Stability stability,
                    std::string output_file_name,
                    const hemelb::lb::GlobalLatticeData &iGlobalLatticeData,
                    const hemelb::lb::LocalLatticeData &iLocalLatticeData);
        void
        WriteConfigParallel(hemelb::lb::Stability stability,
                            std::string output_file_name,
                            const hemelb::lb::GlobalLatticeData &iGlobalLatticeData,
                            const hemelb::lb::LocalLatticeData &iLocalLatticeData);

        double GetMinPhysicalPressure();
        double GetMaxPhysicalPressure();
        double GetMinPhysicalVelocity();
        double GetMaxPhysicalVelocity();
        double GetMinPhysicalStress();
        double GetMaxPhysicalStress();

        void InitMinMaxValues(void);

        double GetAverageInletVelocity(int iInletNumber);
        double GetPeakInletVelocity(int iInletNumber);

        void ReadVisParameters();

        void CalculateMouseFlowField(hemelb::vis::ColPixel *col_pixel_p,
                                     double &mouse_pressure,
                                     double &mouse_stress,
                                     double density_threshold_min,
                                     double density_threshold_minmax_inv,
                                     double stress_threshold_max_inv);

        const hemelb::lb::LbmParameters *GetLbmParams();

      private:
        void CalculateBC(double f[],
                         hemelb::lb::SiteType iSiteType,
                         unsigned int iBoundaryId,
                         double *density,
                         double *vx,
                         double *vy,
                         double *vz,
                         double f_neq[]);

        void InitCollisions();

        //  static void ReadBlock();

        void
        ReadConfig(hemelb::lb::GlobalLatticeData &bGlobalLatticeData);

        void ReadParameters();

        void allocateInlets(int nInlets);
        void allocateOutlets(int nOutlets);

        void handleIOError(int iError);

        double ConvertPressureGradToLatticeUnits(double pressure_grad) const;
        double ConvertPressureGradToPhysicalUnits(double pressure_grad) const;

        hemelb::lb::collisions::MidFluidCollision* mMidFluidCollision;
        hemelb::lb::collisions::WallCollision* mWallCollision;
        hemelb::lb::collisions::InletOutletCollision* mInletCollision;
        hemelb::lb::collisions::InletOutletCollision* mOutletCollision;
        hemelb::lb::collisions::InletOutletWallCollision* mInletWallCollision;
        hemelb::lb::collisions::InletOutletWallCollision* mOutletWallCollision;

        //TODO Get rid of this hack
        hemelb::lb::collisions::Collision* GetCollision(int i);

        double *inlet_density_avg, *inlet_density_amp;
        double *outlet_density_avg, *outlet_density_amp;
        double *inlet_density_phs, *outlet_density_phs;
        int site_min_x, site_min_y, site_min_z;
        int site_max_x, site_max_y, site_max_z;
        int is_inlet_normal_available;
        double* inlet_density, *outlet_density;
        hemelb::lb::collisions::MinsAndMaxes mMinsAndMaxes;
        double *inlet_normal;
        long int *inlet_count;
        double voxel_size;

        double mFileReadTime;

        hemelb::lb::LbmParameters mParams;
        const hemelb::topology::NetworkTopology * mNetTopology;
        hemelb::SimConfig *mSimConfig;

        double *average_inlet_velocity;
        double *peak_inlet_velocity;

        int * receivedFTranslator;
    };
  }
}
#endif // HEMELB_LB_H
