#ifndef HEMELB_LB_H
#define HEMELB_LB_H

#include "net.h"
#include "topology/NetworkTopology.h"
#include "lb/collisions/Collisions.h"
#include "vis/ColPixel.h"
#include "SimConfig.h"

class LBM
{
  public:
    int total_fluid_sites;
    int inlets, outlets;
    int steering_session_id;
    int period;

    double lbmConvertPressureToLatticeUnits(double pressure) const;
    double lbmConvertPressureToPhysicalUnits(double pressure) const;
    double lbmConvertVelocityToLatticeUnits(double velocity) const;
    double lbmConvertStressToLatticeUnits(double stress) const;
    double lbmConvertStressToPhysicalUnits(double stress) const;
    double lbmConvertVelocityToPhysicalUnits(double velocity) const;

    LBM(hemelb::SimConfig *iSimulationConfig,
        const hemelb::topology::NetworkTopology * iNetTop,
        hemelb::lb::GlobalLatticeData &bGlobLatDat,
        int iSteeringSessionId,
        int iPeriod,
        double iVoxelSize,
        Net *net);
    void lbmRestart(hemelb::lb::LocalLatticeData &iLocalLatDat);
    ~LBM();

    int IsUnstable(hemelb::lb::LocalLatticeData &iLocalLatDat);

    hemelb::lb::Stability lbmCycle(int perform_rt,
                                   Net *net,
                                   hemelb::lb::LocalLatticeData &bLocallatDat,
                                   double &bLbTime,
                                   double &bMPISendTime,
                                   double &bMPIWaitTime);
    void lbmCalculateFlowFieldValues();
    void RecalculateTauViscosityOmega();
    void lbmUpdateBoundaryDensities(int cycle_id, int time_step);
    void lbmUpdateInletVelocities(int time_step,
                                  hemelb::lb::LocalLatticeData &iLocalLatDat,
                                  Net *net);

    void lbmSetInitialConditions(hemelb::lb::LocalLatticeData &bLocalLatDat);

    void
    lbmWriteConfig(hemelb::lb::Stability stability,
                   std::string output_file_name,
                   const hemelb::lb::GlobalLatticeData &iGlobalLatticeData,
                   const hemelb::lb::LocalLatticeData &iLocalLatticeData);
    void
        lbmWriteConfigParallel(hemelb::lb::Stability stability,
                               std::string output_file_name,
                               const hemelb::lb::GlobalLatticeData &iGlobalLatticeData,
                               const hemelb::lb::LocalLatticeData &iLocalLatticeData);

    double GetMinPhysicalPressure();
    double GetMaxPhysicalPressure();
    double GetMinPhysicalVelocity();
    double GetMaxPhysicalVelocity();
    double GetMinPhysicalStress();
    double GetMaxPhysicalStress();

    void lbmInitMinMaxValues(void);

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
    void lbmCalculateBC(double f[],
                        hemelb::lb::SiteType iSiteType,
                        unsigned int iBoundaryId,
                        double *density,
                        double *vx,
                        double *vy,
                        double *vz,
                        double f_neq[]);

    void lbmInitCollisions();

    //  static void ReadBlock();

    void
    lbmReadConfig(Net *net, hemelb::lb::GlobalLatticeData &bGlobalLatticeData);

    void lbmReadParameters();

    void allocateInlets(int nInlets);
    void allocateOutlets(int nOutlets);

    void handleIOError(int iError);

    double lbmConvertPressureGradToLatticeUnits(double pressure_grad) const;
    double lbmConvertPressureGradToPhysicalUnits(double pressure_grad) const;

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
    double *lbm_inlet_normal;
    long int *lbm_inlet_count;
    double voxel_size;

    hemelb::lb::LbmParameters mParams;
    const hemelb::topology::NetworkTopology * mNetTopology;
    hemelb::SimConfig *mSimConfig;

    double *lbm_average_inlet_velocity;
    double *lbm_peak_inlet_velocity;

};

#endif // HEMELB_LB_H
