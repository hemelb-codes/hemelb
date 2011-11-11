#include <vector>
#include <cassert>

#include "constants.h"
#include "vis/streaklineDrawer/NeighbouringProcessor.h"

namespace hemelb
{
  namespace vis
  {
    namespace streaklinedrawer
    {
      SendableParticle::SendableParticle() :
        mX(0.0F), mY(0.0F), mZ(0.0F), mVel(0.0F), mInletID(-1)
      {
      }

      SendableParticle::SendableParticle(const Particle& iParticle)
      {
        mX = iParticle.x;
        mY = iParticle.y;
        mZ = iParticle.z;

        mVel = iParticle.vel;

        mInletID = iParticle.inletID;
      }

      Particle SendableParticle::GetParticle()
      {
        return (Particle(mX, mY, mZ, mInletID));
      }

      NeighbouringProcessor::NeighbouringProcessor(proc_t iID) :
        mID(iID)
      {
      }

      void NeighbouringProcessor::AddParticleToSend(const Particle& iParticle)
      {
        mParticlesToSend.push_back(SendableParticle(iParticle));
      }

      bool NeighbouringProcessor::ParticlesToBeRetrieved()
      {
        return (mParticlesToReceive.size() > 0);
      }

      const Particle& NeighbouringProcessor::RetrieveNextReceivedParticle()
      {
        const Particle& lParticle(mParticlesToReceive.back().GetParticle());
        mParticlesToReceive.pop_back();

        return lParticle;
      }

      void NeighbouringProcessor::PrepareToReceiveParticles()
      {
        MPI_Irecv(&mNumberOfParticlesToReceive, 1, //Count
                  MpiDataType<site_t> (), //Type
                  mID, //Destination
                  30, //Tag
                  MPI_COMM_WORLD, //Comm
                  &mReceiveRequest); //Request
      }

      void NeighbouringProcessor::PrepareToSendParticles()
      {
        site_t lParticlesToSend = mParticlesToSend.size();
        MPI_Isend(&lParticlesToSend, 1, // Count
                  MpiDataType<site_t> (), //Type
                  mID, //Destination
                  30, //Tag
                  MPI_COMM_WORLD, // Com
                  &mSendRequest); // Request
      }

      void NeighbouringProcessor::WaitForPreparationToReceiveParticles()
      {
        MPI_Wait(&mReceiveRequest, MPI_STATUS_IGNORE);
      }

      void NeighbouringProcessor::SendParticles()
      {
        if (mParticlesToSend.size() > 0)
        {
          MPI_Isend(&mParticlesToSend[0], (int) mParticlesToSend.size(), //Count
                    MpiDataType<SendableParticle> (), //Type
                    mID, //Destination
                    40, //Tag
                    MPI_COMM_WORLD, //Comm
                    &mSendRequest); //Request
        }
      }

      void NeighbouringProcessor::WaitForParticlesToBeSent()
      {
        if (mParticlesToSend.size() > 0)
        {
          MPI_Wait(&mSendRequest, MPI_STATUS_IGNORE);
        }
        //TODO: need to clear them?
      }

      void NeighbouringProcessor::ReceiveParticles()
      {

        mParticlesToReceive.resize(mNumberOfParticlesToReceive);

        if (mNumberOfParticlesToReceive > 0)
        {
          MPI_Irecv(&mParticlesToReceive[0], (int) mNumberOfParticlesToReceive, //Count
                    MpiDataType<SendableParticle> (), //Type
                    mID, //Source
                    40, //Tag
                    MPI_COMM_WORLD, //Comm
                    &mReceiveRequest); //Request
        }
      }

      void NeighbouringProcessor::WaitForParticlesToBeReceived()
      {
        if (mParticlesToReceive.size() > 0)
        {
          MPI_Wait(&mReceiveRequest, MPI_STATUS_IGNORE);
        }
      }
    }
  }

  template<>
  MPI_Datatype MpiDataTypeTraits<hemelb::vis::streaklinedrawer::SendableParticle>::RegisterMpiDataType()
  {
    int col_pixel_count = 6;
    int col_pixel_blocklengths[6] = { 1, 1, 1, 1, 1, 1 };

    MPI_Datatype col_pixel_types[6] = { MPI_FLOAT,
                                        MPI_FLOAT,
                                        MPI_FLOAT,
                                        MPI_FLOAT,
                                        MPI_UNSIGNED,
                                        MPI_UB };

    MPI_Aint col_pixel_disps[6];

    col_pixel_disps[0] = 0;

    for (int i = 1; i < col_pixel_count; i++)
    {
      if (col_pixel_types[i - 1] == MPI_FLOAT)
      {
        col_pixel_disps[i] = col_pixel_disps[i - 1] + (sizeof(float)
            * col_pixel_blocklengths[i - 1]);
      }
      else if (col_pixel_types[i - 1] == MPI_INT)
      {
        col_pixel_disps[i] = col_pixel_disps[i - 1] + (sizeof(int) * col_pixel_blocklengths[i - 1]);
      }
      else if (col_pixel_types[i - 1] == MPI_UNSIGNED)
      {
        col_pixel_disps[i] = col_pixel_disps[i - 1] + (sizeof(unsigned) * col_pixel_blocklengths[i
            - 1]);
      }

    }
    MPI_Datatype type;
    MPI_Type_struct(col_pixel_count,
                    col_pixel_blocklengths,
                    col_pixel_disps,
                    col_pixel_types,
                    &type);
    MPI_Type_commit(&type);
    return type;
  }
}
