// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#include "colloids/ParticleSet.h"
#include "colloids/BodyForces.h"
#include "colloids/BoundaryConditions.h"
#include <algorithm>
#include "log/Logger.h"
#include "io/writers/xdr/XdrMemWriter.h"
#include "io/formats/formats.h"
#include "io/formats/colloids.h"

namespace hemelb
{
  namespace colloids
  {
    ParticleSet::ParticleSet(const geometry::LatticeData& latDatLBM,
                             io::xml::XmlAbstractionLayer& xml,
                             lb::MacroscopicPropertyCache& propertyCache,
                             std::vector<proc_t>& neighbourProcessors,
                             const std::string& outputPath) :
        localRank(topology::NetworkTopology::Instance()->GetLocalRank()), latDatLBM(latDatLBM), propertyCache(propertyCache), path(outputPath)
    {
      // add an element into scanMap for each neighbour rank with zero for both counts
      // sorting the list of neighbours allows the position in the map to be predicted
      // & giving the correct position in the map makes insertion significantly faster
      // the local rank is added last, because its position cannot be easily predicted
      std::sort(neighbourProcessors.begin(), neighbourProcessors.end());
      for (std::vector<proc_t>::const_iterator iter = neighbourProcessors.begin(); iter != neighbourProcessors.end();
          iter++)
        scanMap.insert(scanMap.end(), scanMapContentType(*iter, scanMapElementType(0, 0)));
      scanMap.insert(scanMapContentType(localRank, scanMapElementType(0, 0)));

      // assume we are at the <Particles> node
      bool found = xml.MoveToChild("subgridParticle");
      if (found)
        propertyCache.velocityCache.SetRefreshFlag();
      while (found)
      {
        // create the particle object from the settings in the config file
        Particle nextParticle(latDatLBM, xml);
        // check the particle is valid, i.e. in fluid, and is locally owned
        if (nextParticle.IsValid() && nextParticle.GetOwnerRank() == localRank)
        {
          // add the particle to the list of known particles ...
          particles.push_back(nextParticle);
          // ... and keep the count of local particles up-to-date
          scanMap[localRank].first++;
        }
        found = xml.NextSibling("subgridParticle");
      }
    }
    ;

    ParticleSet::~ParticleSet()
    {
      particles.clear();
    }

    const void ParticleSet::OutputInformation(const LatticeTime timestep)
    {
      const unsigned int maxSize = io::formats::colloids::RecordLength * particles.size()
          + io::formats::colloids::HeaderLength + io::formats::colloids::MagicLength;
      if (buffer.size() < maxSize)
      {
        buffer.resize(maxSize);
      }

      io::writers::xdr::XdrMemWriter writer(&buffer.front(), maxSize);

      for (std::vector<Particle>::iterator iter = particles.begin(); iter != particles.end(); iter++)
      {
        Particle& particle = *iter;
        if (particle.GetOwnerRank() == localRank)
        {
          particle.OutputInformation();
          particle.WriteToStream(timestep, * ((io::writers::Writer*) &writer));
        }
      }
      const unsigned int count = writer.getCurrentStreamPosition();

      // use MPI-IO to write the buffer to colloid output file
      MPI_File fh;
      MPI_File_open(MPI_COMM_WORLD,
                    const_cast<char*>(path.c_str()),
                    MPI_MODE_APPEND | MPI_MODE_WRONLY | MPI_MODE_CREATE,
                    MPI_INFO_NULL,
                    &fh);

      // work-around: the shared file pointer may not be set correctly
      //              on all ranks immediately after opening the file,
      //              or indeed after any shared/collective operation.
      //              In particular this happens on martensite at EPCC
      //              using MPICH2-1.4.1p1 on 64bit linux 2.6.18 SL5.0
      //              Issuing "sync;barrier;sync" enforces consistency

      MPI_File_sync(fh);
      MPI_Barrier (MPI_COMM_WORLD);
      MPI_File_sync(fh);

      MPI_Offset offsetEOF;
      MPI_File_get_position_shared(fh, &offsetEOF);

      MPI_File_sync(fh);
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_File_sync(fh);

      MPI_Offset dispStartOfHeader;
      MPI_File_get_byte_offset(fh, offsetEOF, &dispStartOfHeader);

      unsigned int sizeOfHeader = io::formats::colloids::HeaderLength;
      if (dispStartOfHeader == 0)
      {
        // the header starts at the begining of a new file, so we
        // write the magic numbers that identify the type of file
        sizeOfHeader += io::formats::colloids::MagicLength;
      }

      log::Logger::Log<log::Debug, log::OnePerCore>("dispStartOfHeader: %i (from offsetEOF: %i)\n",
                                                    dispStartOfHeader,
                                                    offsetEOF);

      MPI_File_set_view(fh, dispStartOfHeader + sizeOfHeader, MPI_CHAR, MPI_CHAR, "native\0", MPI_INFO_NULL);

      log::Logger::Log<log::Debug, log::OnePerCore>("SetView for data - disp: %i\n", dispStartOfHeader + sizeOfHeader);

      MPI_File_sync(fh);
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_File_sync(fh);

      // collective write: the effect is as though all writes are done
      // in serialised order, i.e. as if rank 0 writes first, followed
      // by rank 1, and so on, until all ranks have written their data
      MPI_File_write_ordered(fh, &buffer.front(), count, MPI_CHAR, MPI_STATUS_IGNORE);

      MPI_File_sync(fh);
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_File_sync(fh);

      // the collective ordered write modifies the shared file pointer
      // it should point to the byte following the highest rank's data
      // (should be true for all ranks but) we only need it for rank 0
      MPI_File_get_position_shared(fh, &offsetEOF);

      MPI_File_sync(fh);
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_File_sync(fh);

      // only rank 0 uses this view but this is a collective operation
      MPI_File_set_view(fh, dispStartOfHeader, MPI_CHAR, MPI_CHAR, "native\0", MPI_INFO_NULL);

      log::Logger::Log<log::Debug, log::OnePerCore>("dispStartOfHeader: %i (new offsetEOF: %i)\n",
                                                    dispStartOfHeader,
                                                    offsetEOF);

      MPI_File_sync(fh);
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_File_sync(fh);

      if (localRank == 0)
      {
        if (dispStartOfHeader == 0)
        {
          writer << (uint32_t) io::formats::HemeLbMagicNumber;
          writer << (uint32_t) io::formats::colloids::MagicNumber;
          writer << (uint32_t) io::formats::colloids::VersionNumber;
        }
        writer << (uint32_t) io::formats::colloids::HeaderLength;
        writer << (uint32_t) io::formats::colloids::RecordLength;
        writer << (uint64_t) offsetEOF;
        writer << (uint64_t) timestep;
        MPI_File_write(fh, &buffer[count], sizeOfHeader, MPI_CHAR, MPI_STATUS_IGNORE);
      }

      MPI_File_close(&fh);

      for (scanMapConstIterType iterMap = scanMap.begin(); iterMap != scanMap.end(); iterMap++)
      {
        const proc_t& neighbourRank = iterMap->first;
        const unsigned int& numberOfParticles = iterMap->second.first;
        const unsigned int& numberOfVelocities = iterMap->second.second;
        log::Logger::Log<log::Debug, log::OnePerCore>("ScanMap[%i] = {%i, %i}\n",
                                                      neighbourRank,
                                                      numberOfParticles,
                                                      numberOfVelocities);
      }
    }

    const void ParticleSet::UpdatePositions()
    {
      if (log::Logger::ShouldDisplay<log::Debug>())
        log::Logger::Log<log::Debug, log::OnePerCore>("In colloids::ParticleSet::UpdatePositions #particles == %i ...\n",
                                                      topology::NetworkTopology::Instance()->GetLocalRank(),
                                                      particles.size());

      // only update the position for particles that are locally owned because
      // only the owner has velocity contributions from all neighbouring ranks
      for (std::vector<Particle>::iterator iter = particles.begin(); iter != particles.end(); iter++)
      {
        Particle& particle = *iter;
        if (particle.GetOwnerRank() == localRank)
          particle.UpdatePosition(latDatLBM);
      }
    }

    const void ParticleSet::CalculateBodyForces()
    {
      for (std::vector<Particle>::iterator iter = particles.begin(); iter != particles.end(); iter++)
      {
        Particle& particle = *iter;
        if (particle.GetOwnerRank() == localRank)
          particle.CalculateBodyForces();
      }
    }

    const void ParticleSet::ApplyBoundaryConditions(const LatticeTime currentTimestep)
    {
      for (std::vector<Particle>::iterator iter = particles.begin(); iter != particles.end(); iter++)
      {
        Particle& particle = *iter;
        if (particle.GetOwnerRank() == localRank)
        {
          BoundaryConditions::DoSomeThingsToParticle(currentTimestep, particle);
          if (particle.IsReadyToBeDeleted())
            log::Logger::Log<log::Trace, log::OnePerCore>("In ParticleSet::ApplyBoundaryConditions - timestep: %lu, particleId: %lu, IsReadyToBeDeleted: %s, markedForDeletion: %lu, lastCheckpoint: %lu\n",
                                                          currentTimestep,
                                                          particle.GetParticleId(),
                                                          particle.IsReadyToBeDeleted() ?
                                                            "YES" :
                                                            "NO",
                                                          particle.GetDeletionMarker(),
                                                          particle.GetLastCheckpointTimestep());
        }
      }

      // shuffle (or partition) the particles in our vector containing all particles
      // so the first partition contains all the local particles that should be kept
      // and the other contains all the deletable-local plus the non-local particles
      std::vector<Particle>::iterator bound =
          std::partition(particles.begin(),
                         particles.begin() + scanMap[localRank].first,
                         std::not1(std::mem_fun_ref(&Particle::IsReadyToBeDeleted)));

      if (scanMap[localRank].first > (bound - particles.begin()))
        log::Logger::Log<log::Debug, log::OnePerCore>("In ParticleSet::ApplyBoundaryConditions - timestep: %lu, scanMap[localRank].first: %lu, bound-particles.begin(): %lu\n",
                                                      currentTimestep,
                                                      scanMap[localRank].first,
                                                      bound - particles.begin());

      // the partitioning above may invalidate the scanMap used by the communication
      // the next communication function called is CommunicatePositions, which needs
      // the number of local particles to be correct - i.e. scanMap[localRank].first
      // - the rest of scanMap will be re-built by the CommunicatePositions function
      scanMap[localRank].first = bound - particles.begin();
    }

    const void ParticleSet::CalculateFeedbackForces()
    {
      BodyForces::ClearBodyForcesForAllSiteIds();
      for (std::vector<Particle>::const_iterator iter = particles.begin(); iter != particles.end(); iter++)
      {
        const Particle& particle = *iter;
        if (particle.GetOwnerRank() == localRank)
          particle.CalculateFeedbackForces(latDatLBM);
      }
    }

    const void ParticleSet::InterpolateFluidVelocity()
    {
      for (std::vector<Particle>::iterator iter = particles.begin(); iter != particles.end(); iter++)
      {
        Particle& particle = *iter;
        particle.InterpolateFluidVelocity(latDatLBM, propertyCache);
      }
      propertyCache.velocityCache.SetRefreshFlag();
    }

    const void ParticleSet::CommunicateParticlePositions()
    {
      /** CommunicateParticlePositions
       *    For each neighbour rank p
       *    - MPI_IRECV( number_of_remote_particles )
       *    - MPI_IRECV( list_of_remote_particles )
       *    - MPI_ISEND( number_of_local_particles )
       *    - MPI_ISEND( list_of_local_particles )
       *    MPI_WAITALL()
       *
       *  The global position of each particle is updated by the ownerRank process.
       *  The ownerRank for each particle is verified when its position is updated.
       *  Some (previously locally owned) particles may no longer be locally owned.
       *  
       */

      unsigned int& numberOfParticlesToSend = scanMap[localRank].first;
      if (scanMap.size() < 2)
      {
        return;
      }

      for (scanMapIterType iterMap = scanMap.begin(); iterMap != scanMap.end(); iterMap++)
      {
        const proc_t& neighbourRank = iterMap->first;
        if (neighbourRank != localRank)
        {
          unsigned int& numberOfParticlesToRecv = iterMap->second.first;
          net.RequestSendR(numberOfParticlesToSend, neighbourRank);
          net.RequestReceiveR(numberOfParticlesToRecv, neighbourRank);
        }
      }
      net.Dispatch();

      unsigned int numberOfParticles = 0;
      for (scanMapConstIterType iterMap = scanMap.begin(); iterMap != scanMap.end(); iterMap++)
        numberOfParticles += iterMap->second.first;
      particles.resize(numberOfParticles);

      std::vector<Particle>::iterator iterSendBegin = particles.begin();
      std::vector<Particle>::iterator iterRecvBegin = particles.begin() + numberOfParticlesToSend;
      for (scanMapConstIterType iterMap = scanMap.begin(); iterMap != scanMap.end(); iterMap++)
      {
        const proc_t& neighbourRank = iterMap->first;
        if (neighbourRank != localRank)
        {
          const unsigned int& numberOfParticlesToRecv = iterMap->second.first;
          net.RequestSend(& ((PersistedParticle&) *iterSendBegin), numberOfParticlesToSend, neighbourRank);
          net.RequestReceive(& ((PersistedParticle&) * (iterRecvBegin)), numberOfParticlesToRecv, neighbourRank);
          iterRecvBegin += numberOfParticlesToRecv;
        }
      }
      net.Dispatch();

      // remove particles owned by unknown ranks
      std::vector<Particle>::iterator newEndOfParticles =
          std::partition(particles.begin(),
                         particles.end(),
                         std::bind2nd(std::mem_fun_ref(&Particle::IsOwnerRankKnown), scanMap));
      particles.erase(newEndOfParticles, particles.end());

      // sort the particles - local first, then in order of increasing owner rank
      std::sort(particles.begin(), particles.end());

      // re-build the scanMap
      for (scanMapIterType iterMap = scanMap.begin(); iterMap != scanMap.end(); iterMap++)
        iterMap->second.first = 0;
      for (std::vector<Particle>::const_iterator iterParticles = particles.begin(); iterParticles != particles.end();
          iterParticles++)
        scanMap[iterParticles->GetOwnerRank()].first++;
    }

    const void ParticleSet::CommunicateFluidVelocities()
    {
      /** CommunicateFluidVelocities
       *    For each neighbour rank p
       *    - MPI_IRECV( number_of_incoming_velocities )
       *    - MPI_IRECV( list_of_incoming_velocities )
       *    - MPI_ISEND( number_of_outgoing_velocities )
       *    - MPI_ISEND( list_of_outgoing_velocities )
       *    MPI_WAITALL()
       */

      if (scanMap.size() < 2)
      {
        return;
      }

      // exchange counts
      for (scanMapIterType iterMap = scanMap.begin(); iterMap != scanMap.end(); iterMap++)
      {
        const proc_t& neighbourRank = iterMap->first;
        if (neighbourRank != localRank)
        {
          unsigned int& numberOfVelocitiesToSend = iterMap->second.first;
          unsigned int& numberOfVelocitiesToRecv = iterMap->second.second;
          net.RequestSendR(numberOfVelocitiesToSend, neighbourRank);
          net.RequestReceiveR(numberOfVelocitiesToRecv, neighbourRank);
        }
      }
      net.Dispatch();

      // sum counts
      unsigned int numberOfIncomingVelocities = 0;
      for (scanMapConstIterType iterMap = scanMap.begin(); iterMap != scanMap.end(); iterMap++)
        numberOfIncomingVelocities += iterMap->second.second;
      velocityBuffer.resize(numberOfIncomingVelocities);

      // exchange velocities
      std::vector<Particle>::iterator iterSendBegin = particles.begin();
      std::vector<std::pair<unsigned long, util::Vector3D<double> > >::iterator iterRecvBegin = velocityBuffer.begin();
      for (scanMapConstIterType iterMap = scanMap.begin(); iterMap != scanMap.end(); iterMap++)
      {
        const proc_t& neighbourRank = iterMap->first;
        if (neighbourRank != localRank)
        {
          const unsigned int& numberOfVelocitiesToSend = iterMap->second.first;
          const unsigned int& numberOfVelocitiesToRecv = iterMap->second.second;
          net.RequestSend(& ((Particle&) *iterSendBegin), numberOfVelocitiesToSend, neighbourRank);
          net.RequestReceive(& (* (iterRecvBegin)), numberOfVelocitiesToRecv, neighbourRank);
          iterRecvBegin += numberOfVelocitiesToRecv;
        }
      }
      net.Dispatch();

      // sum velocities
      velocityMap.clear();
      for (std::vector<std::pair<unsigned long, util::Vector3D<double> > >::const_iterator iterVelocityBuffer =
          velocityBuffer.begin(); iterVelocityBuffer != velocityBuffer.end(); iterVelocityBuffer++)
      {
        const unsigned long& particleId = iterVelocityBuffer->first;
        const util::Vector3D<double>& partialVelocity = iterVelocityBuffer->second;
        velocityMap[particleId] += partialVelocity;
      }

      // update local particles
      for (std::vector<Particle>::iterator iter = particles.begin(); iter != particles.end(); iter++)
      {
        Particle& particle = *iter;
        if (particle.GetOwnerRank() == localRank)
          particle.AccumulateVelocity(velocityMap[particle.GetParticleId()]);
      }

    }

  }
}
