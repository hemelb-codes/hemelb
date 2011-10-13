#ifndef HEMELB_VIS_STREAKLINEDRAWER_H
#define HEMELB_VIS_STREAKLINEDRAWER_H

#include <vector>
#include <map>

#include "constants.h"
#include "mpiInclude.h"

#include "geometry/LatticeData.h"
#include "topology/NetworkTopology.h"

#include "vis/PixelSetStore.h"
#include "vis/PixelSet.h"
#include "vis/BasicPixel.h"
#include "vis/Screen.h"
#include "vis/StreakPixel.h"
#include "vis/Viewpoint.h"
#include "vis/VisSettings.h"

namespace hemelb
{
  namespace vis
  {
    //TODO Some of this class's members could be combined with the topology-aware classes
    // (and obv be initialised when they are).

    /**
     * Class that controls the drawing of streaklines - lines that trace
     * the path of an imaginary particle were it dropped into the fluid.
     */
    class StreaklineDrawer : public PixelSetStore<PixelSet<StreakPixel> >
    {
      public:
        // Constructor and destructor.
        StreaklineDrawer(geometry::LatticeData* iLatDat,
                         Screen* iScreen,
                         Viewpoint* iViewpoint,
                         VisSettings* iVisSettings);
        ~StreaklineDrawer();

        // Method to reset streakline drawer
        void Restart();

        // Drawing methods.
        void StreakLines(unsigned long time_steps,
                         unsigned long time_steps_per_cycle,
                         geometry::LatticeData* iLatDat);
        PixelSet<StreakPixel>* Render(geometry::LatticeData* iLatDat);

      private:
        // Struct for a particle dropped into the fluid.
        struct Particle
        {
            float x, y, z;
            float vx, vy, vz;
            float vel;
            unsigned int inlet_id;
        };

        // Struct for information about the velocity field at some point.
        struct VelSiteData
        {
            proc_t proc_id;
            site_t counter, site_id;
            float vx, vy, vz;
        };

        // Struct to contain the whole velocity field.
        struct VelocityField
        {
            VelSiteData *vel_site_data;
        };

        // Struct for one processor to hold information about its neighbouring processors.
        typedef std::vector<float> FloatVector;
        struct NeighProc
        {
            proc_t id;
            site_t send_ps, recv_ps;
            site_t send_vs, recv_vs;

            FloatVector p_to_send, p_to_recv;
            float *v_to_send, *v_to_recv;

            site_t *s_to_send, *s_to_recv;
        };

        // Necessary to keep a local store of the number of blocks created, so that we can
        // write a correct constructor. Alternative (TODO) is to use a vector.
        site_t num_blocks;

        // Counter keeps track of the number of VelSiteDatas created
        site_t counter;

        // Variables for tracking the actual numbers of particles, and the maximum number
        //(i.e. the number for which memory has been allocated).
        unsigned int nParticles;
        unsigned int nParticleSeeds;
        unsigned int particles_to_send_max, particles_to_recv_max;

        // Variables for counting the processors involved etc.
        site_t shared_vs;
        proc_t procs;

        // Pointers to the structs.
        VelocityField *velocity_field;
        std::vector<Particle> particleVec;
        std::vector<Particle> particleSeedVec;

        // Arrays for communicating between processors.
        float *v_to_send, *v_to_recv;
        site_t *s_to_send, *s_to_recv;
        proc_t *from_proc_id_to_neigh_proc_index;

        std::vector<NeighProc> mNeighProcs;

        // Require these for inter-processor comms.
        MPI_Request *req;

        // Private functions for the creation / deletion of particles.
        void createSeedParticles();
        void createParticle(float x, float y, float z, float vel, int inlet_id);
        void deleteParticle(unsigned int p_index);

        // Private functions for initialising the velocity field.
        void initializeVelFieldBlock(const geometry::LatticeData* iLatDat,
                                     site_t site_i,
                                     site_t site_j,
                                     site_t site_k,
                                     proc_t proc_id);
        VelSiteData *velSiteDataPointer(geometry::LatticeData* iLatDat,
                                        site_t site_i,
                                        site_t site_j,
                                        site_t site_k);
        void particleVelocity(Particle *particle_p, float v[2][2][2][3], float interp_v[3]);
        void localVelField(int p_index,
                           float v[2][2][2][3],
                           int *is_interior,
                           geometry::LatticeData* iLatDat);

        // Private functions for updating the velocity field and the particles in it.
        void updateVelField(int stage_id, geometry::LatticeData* iLatDat);
        void updateParticles();

        // Private functions for inter-proc communication.
        void communicateSiteIds();
        void communicateVelocities(geometry::LatticeData* iLatDat);
        void communicateParticles(geometry::LatticeData* iLatDat);

        Screen* mScreen;
        Viewpoint* mViewpoint;
        VisSettings* mVisSettings;
    };

  }
}

#endif // HEMELB_VIS_STREAKLINEDRAWER_H
