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

      NeighbouringProcessor::NeighbouringProcessor(proc_t iID) :
          mID(iID)
      {
      }

      void NeighbouringProcessor::AddParticleToSend(const Particle& iParticle)
      {
        mParticlesToSend.push_back(iParticle);
      }

      bool NeighbouringProcessor::ParticlesToBeRetrieved()
      {
        return (mParticlesToReceive.size() > 0);
      }

      const Particle& NeighbouringProcessor::RetrieveNextReceivedParticle()
      {
        const Particle& lParticle = mParticlesToReceive.back();
        mParticlesToReceive.pop_back();
        return lParticle;
      }

      void NeighbouringProcessor::PrepareToReceiveParticles()
      {
        MPI_Irecv(&mNumberOfParticlesToReceive, 1, //Count
                  MpiDataType<site_t>(), //Type
                  mID, //Destination
                  30, //Tag
                  MPI_COMM_WORLD, //Comm
                  &mReceiveRequest); //Request
      }

      void NeighbouringProcessor::PrepareToSendParticles()
      {
        site_t particlesToSend = mParticlesToSend.size();
        MPI_Isend(&particlesToSend, 1, // Count
                  MpiDataType<site_t>(), //Type
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
          MPI_Isend(&mParticlesToSend[0], mParticlesToSend.size(), //Count
                    MpiDataType<Particle>(), //Type
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
        mParticlesToSend.clear();
      }

      void NeighbouringProcessor::ReceiveParticles()
      {

        mParticlesToReceive.resize(mNumberOfParticlesToReceive);

        if (mNumberOfParticlesToReceive > 0)
        {
          MPI_Irecv(&mParticlesToReceive[0], mNumberOfParticlesToReceive, //Count
                    MpiDataType<Particle>(), //Type
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
  MPI_Datatype MpiDataTypeTraits<hemelb::vis::streaklinedrawer::Particle>::RegisterMpiDataType()
  {
    const int elementCount = 7;
    int elementBlockLengths[elementCount] = { 1, 1, 1, 1, 1, 1, 1 };

    MPI_Datatype elementTypes[elementCount] = { MPI_LB,
                                                    MPI_FLOAT,
                                                    MPI_FLOAT,
                                                    MPI_FLOAT,
                                                    MPI_FLOAT,
                                                    MPI_UNSIGNED,
                                                    MPI_UB };

    MPI_Aint elementDisplacements[elementCount];

    vis::streaklinedrawer::Particle particle[2];

    MPI_Address(&particle[0], &elementDisplacements[0]);
    MPI_Address(&particle[0].x, &elementDisplacements[1]);
    MPI_Address(&particle[0].y, &elementDisplacements[2]);
    MPI_Address(&particle[0].z, &elementDisplacements[3]);
    MPI_Address(&particle[0].vel, &elementDisplacements[4]);
    MPI_Address(&particle[0].inletID, &elementDisplacements[5]);
    MPI_Address(&particle[1], &elementDisplacements[6]);
     for (int element = elementCount - 1; element >= 0; element--) {
       elementDisplacements[element] -= elementDisplacements[0];
     }


    MPI_Datatype type;
    MPI_Type_struct(elementCount, elementBlockLengths, elementDisplacements, elementTypes, &type);
    MPI_Type_commit(&type);
    return type;
  }
}
