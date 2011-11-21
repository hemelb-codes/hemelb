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

      NeighbouringProcessor::NeighbouringProcessor()
      {
      }

      NeighbouringProcessor::NeighbouringProcessor(proc_t neighbourRankIn) :
        neighbourRank(neighbourRankIn)
      {
      }

      void NeighbouringProcessor::AddParticleToSend(const Particle& particle)
      {
        particlesToSend.push_back(particle);
      }

      bool NeighbouringProcessor::ParticlesToBeRetrieved()
      {
        return (particlesToReceive.size() > 0);
      }

      const Particle& NeighbouringProcessor::PopNextReceivedParticle()
      {
        const Particle& lParticle = particlesToReceive.back();
        particlesToReceive.pop_back();
        return lParticle;
      }

      void NeighbouringProcessor::PrepareToReceiveParticles()
      {
        MPI_Irecv(&numberOfParticlesToReceive, 1, //Count
                  MpiDataType<> (numberOfParticlesToReceive), //Type
                  neighbourRank, //Destination
                  30, //Tag
                  MPI_COMM_WORLD, //Comm
                  &receiveRequest); //Request
      }

      void NeighbouringProcessor::PrepareToSendParticles()
      {
        numberOfParticlesToSend = particlesToSend.size();
        MPI_Isend(&numberOfParticlesToSend, 1, // Count
                  MpiDataType<> (numberOfParticlesToSend), //Type
                  neighbourRank, //Destination
                  30, //Tag
                  MPI_COMM_WORLD, // Com
                  &sendRequest); // Request
      }

      void NeighbouringProcessor::WaitForPreparationToReceiveParticles()
      {
        MPI_Wait(&receiveRequest, MPI_STATUS_IGNORE);
      }

      void NeighbouringProcessor::SendParticles()
      {
        if (particlesToSend.size() > 0)
        {
          MPI_Isend(&particlesToSend[0], (int) particlesToSend.size(), //Count
                    MpiDataType<Particle> (), //Type
                    neighbourRank, //Destination
                    40, //Tag
                    MPI_COMM_WORLD, //Comm
                    &sendRequest); //Request
        }
      }

      void NeighbouringProcessor::WaitForParticlesToBeSent()
      {
        if (particlesToSend.size() > 0)
        {
          MPI_Wait(&sendRequest, MPI_STATUS_IGNORE);
        }
        particlesToSend.clear();
      }

      void NeighbouringProcessor::ReceiveParticles()
      {
        particlesToReceive.resize(numberOfParticlesToReceive);

        if (numberOfParticlesToReceive > 0)
        {
          MPI_Irecv(&particlesToReceive[0], (int) numberOfParticlesToReceive, //Count
                    MpiDataType<Particle> (), //Type
                    neighbourRank, //Source
                    40, //Tag
                    MPI_COMM_WORLD, //Comm
                    &receiveRequest); //Request
        }
      }

      void NeighbouringProcessor::WaitForParticlesToBeReceived()
      {
        if (particlesToReceive.size() > 0)
        {
          MPI_Wait(&receiveRequest, MPI_STATUS_IGNORE);
        }
      }

      void NeighbouringProcessor::AddSiteToRequestVelocityDataFor(site_t siteI,
                                                                  site_t siteJ,
                                                                  site_t siteK)
      {
        siteCoordsToSend.push_back(util::Vector3D<site_t>(siteI, siteJ, siteK));
      }

      void NeighbouringProcessor::ExchangeSiteIds()
      {
        sentSitesCount = siteCoordsToSend.size();
        MPI_Irecv(&receivedSitesCount,
                  1,
                  MpiDataType(receivedSitesCount),
                  neighbourRank,
                  30,
                  MPI_COMM_WORLD,
                  &siteIdReceiveRequest);
        MPI_Isend(&sentSitesCount,
                  1,
                  MpiDataType(sentSitesCount),
                  neighbourRank,
                  30,
                  MPI_COMM_WORLD,
                  &siteIdSendRequest);

        MPI_Wait(&siteIdReceiveRequest, MPI_STATUS_IGNORE);
        MPI_Wait(&siteIdSendRequest, MPI_STATUS_IGNORE);

        if (receivedSitesCount > 0)
        {
          siteCoordsToReceive.resize(receivedSitesCount);
          sentVelocityField.resize(receivedSitesCount);

          MPI_Irecv(&siteCoordsToReceive[0],
                    (int) receivedSitesCount,
                    MpiDataType(siteCoordsToReceive[0]),
                    neighbourRank,
                    40,
                    MPI_COMM_WORLD,
                    &siteIdReceiveRequest);
        }

        if (sentSitesCount > 0)
        {
          receivedVelocityField.resize(sentSitesCount);
          MPI_Isend(&siteCoordsToSend[0],
                    (int) sentSitesCount,
                    MpiDataType(siteCoordsToSend[0]),
                    neighbourRank,
                    40,
                    MPI_COMM_WORLD,
                    &siteIdSendRequest);
        }
      }

      void NeighbouringProcessor::WaitForSiteIdExchange()
      {
        if (receivedSitesCount > 0)
        {
          MPI_Wait(&siteIdReceiveRequest, MPI_STATUS_IGNORE);
        }

        if (sentSitesCount > 0)
        {
          MPI_Wait(&siteIdSendRequest, MPI_STATUS_IGNORE);
        }
      }

      void NeighbouringProcessor::ExchangeVelocitiesForRequestedSites()
      {
        if (sentSitesCount > 0)
        {
          receivedVelocityField.resize(sentSitesCount);

          MPI_Irecv(&receivedVelocityField[0],
                    (int) sentSitesCount,
                    MpiDataType(receivedVelocityField[0]),
                    neighbourRank,
                    30,
                    MPI_COMM_WORLD,
                    &velocityReceiveRequest);
        }

        if (receivedSitesCount > 0)
        {
          sentVelocityField.resize(receivedSitesCount);

          MPI_Isend(&sentVelocityField[0],
                    (int) receivedSitesCount,
                    MpiDataType(sentVelocityField[0]),
                    neighbourRank,
                    30,
                    MPI_COMM_WORLD,
                    &velocitySendRequest);
        }
      }

      void NeighbouringProcessor::WaitForVelocityExchange()
      {
        if (sentSitesCount > 0)
        {
          MPI_Wait(&velocityReceiveRequest, MPI_STATUS_IGNORE);
        }
        if (receivedSitesCount > 0)
        {
          MPI_Wait(&velocitySendRequest, MPI_STATUS_IGNORE);
        }
      }

      site_t NeighbouringProcessor::GetSitesToReceiveCount()
      {
        return siteCoordsToReceive.size();
      }

      const util::Vector3D<float>& NeighbouringProcessor::GetReceivedVelocity(const size_t receivedIndex) const
      {
        return receivedVelocityField[receivedIndex];
      }

      const util::Vector3D<site_t>& NeighbouringProcessor::GetSiteCoordsBeingReceived(const size_t receivedIndex) const
      {
        return siteCoordsToReceive[receivedIndex];
      }

      void NeighbouringProcessor::SetVelocityFieldToSend(const size_t sendIndex,
                                                         const util::Vector3D<float>& velocityFieldToSend)
      {
        sentVelocityField[sendIndex] = velocityFieldToSend;
      }

      size_t NeighbouringProcessor::GetSendingSiteCount() const
      {
        return siteCoordsToSend.size();
      }

      util::Vector3D<site_t>& NeighbouringProcessor::GetSiteCoorinates(size_t sendIndex)
      {
        return siteCoordsToSend[sendIndex];
      }

      void NeighbouringProcessor::ClearSendingSites()
      {
        siteCoordsToSend.clear();
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
    for (int element = elementCount - 1; element >= 0; element--)
    {
      elementDisplacements[element] -= elementDisplacements[0];
    }

    MPI_Datatype type;
    MPI_Type_struct(elementCount, elementBlockLengths, elementDisplacements, elementTypes, &type);
    MPI_Type_commit(&type);
    return type;
  }
}
