#ifndef HEMELB_VIS_STREAKLINEDRAWER_NEIGHBOURINGPROCESSOR_H
#define HEMELB_VIS_STREAKLINEDRAWER_NEIGHBOURINGPROCESSOR_H

#include <vector>

#include "mpiInclude.h"
#include "vis/streaklineDrawer/Particle.h"

namespace hemelb
{
  namespace vis
  {
    namespace streaklinedrawer
    {
      class SendableParticle
      {
        public:
          SendableParticle();
          SendableParticle(const Particle& iParticle);

          Particle GetParticle();

        private:
          float mX;
          float mY;
          float mZ;

          float mVel;

          unsigned int mInletID;
      };

      class NeighbouringProcessor
      {
        public:
          NeighbouringProcessor(proc_t iID);

          void AddParticleToSend(const Particle& iParticle);

          bool ParticlesToBeRetrieved();

          const Particle& RetrieveNextReceivedParticle();

          void PrepareToReceiveParticles();

          void PrepareToSendParticles();

          void WaitForPreparationToReceiveParticles();

          void SendParticles();

          void WaitForParticlesToBeSent();

          void ReceiveParticles();

          void WaitForParticlesToBeReceived();

          proc_t mID;

          site_t send_vs, recv_vs;

          float *v_to_send, *v_to_recv;
          site_t *s_to_send, *s_to_recv;

        private:
          std::vector<SendableParticle> mParticlesToSend;

          size_t mNumberOfParticlesToReceive;
          std::vector<SendableParticle> mParticlesToReceive;

          MPI_Request mSendRequest;
          MPI_Request mReceiveRequest;
      };
    }
  }
}

#endif // HEMELB_VIS_STREAKLINEDRAWER_NEIGHBOURINGPROCESSOR_H
