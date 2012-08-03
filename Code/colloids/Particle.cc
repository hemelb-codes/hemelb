#include "colloids/Particle.h"
#include "geometry/LatticeData.h"
#include "log/Logger.h"
#include "units.h"

namespace hemelb
{
  // boiler-plate template specialisation for colloids Particle object
  template<>
  MPI_Datatype MpiDataTypeTraits<colloids::Particle>::RegisterMpiDataType()
  {
    return colloids::Particle().CreateMpiDatatypeWithVelocity();
  }

  // boiler-plate template specialisation for colloids PersistedParticle object
  template<>
  MPI_Datatype MpiDataTypeTraits<colloids::PersistedParticle>::RegisterMpiDataType()
  {
    return colloids::Particle().CreateMpiDatatypeWithPosition();
  }

  // boiler-plate template specialisation for incoming particle velocity
  template<>
  MPI_Datatype MpiDataTypeTraits<std::pair<unsigned long, util::Vector3D<double> > >::RegisterMpiDataType()
  {
    // MPI::Get_address specifies non-const pointers
    // so need to create a non-const std::pair object
    std::pair<unsigned long, util::Vector3D<double> > temp;

    // we want a re-usable MPI data type
    // so we need relative displacements
    MPI::Aint baseAddress = MPI::Get_address(&temp);

    // we have chosen to make each block of fields contain a single field
    // so, the number of field blocks is the same as the number of fields
    int numberOfFieldBlocks = 2;

    // and the length of every field block is one
    int lengthOfEachFieldBlock[] = {1, 1};

    // there is no guarantee that the fields will be contiguous, so
    // the displacement of each field must be determined separately
    MPI::Aint displacementOfEachFieldBlock[numberOfFieldBlocks];
    displacementOfEachFieldBlock[0] = MPI::Get_address(&(temp.first)) - baseAddress;
    displacementOfEachFieldBlock[1] = MPI::Get_address(&(temp.second)) - baseAddress;

    // the built-in MPI datatype of each field must match the C++ type
    MPI::Datatype datatypeOfEachFieldBlock[] =
      {MPI::UNSIGNED_LONG, MpiDataType<util::Vector3D<double> >()};

    // create a first draft of the MPI datatype for a Particle
    // the lower bound and displacements of fields are correct
    // but the extent may not include the whole derived object
    // specifically, we aren't sending velocity and bodyForces
    MPI::Datatype pairType = MPI::Datatype::Create_struct(
      numberOfFieldBlocks,
      lengthOfEachFieldBlock,
      displacementOfEachFieldBlock,
      datatypeOfEachFieldBlock);

    // obtain the current lower bound for the MPI datatype
    MPI::Aint lowerBound, extent;
    pairType.Get_extent(lowerBound, extent);

    // we can determine the actual extent of a Particle object
    // by concatenating two of them, using a contiguous vector
    // and finding the difference between their base addresses
    std::vector<std::pair<unsigned long, util::Vector3D<double> > > tempVectorOfPair(2, temp);
    extent = MPI::Get_address(&(tempVectorOfPair[1]))
           - MPI::Get_address(&(tempVectorOfPair[0]));

    // resize the uncommitted first draft MPI datatype
    // with the current lower bound and the new extent
    pairType = pairType.Create_resized(lowerBound, extent);

    // commit the MPI datatype and return it
    pairType.Commit();
    return pairType;
  }

  namespace colloids
  {
    Particle::Particle(const geometry::LatticeData& latDatLBM,
                       io::xml::XmlAbstractionLayer& xml) :
      PersistedParticle(xml)
    {
      /** TODO: this conversion should be done in the xml abstraction layer */
      smallRadius_a0 /= latDatLBM.GetVoxelSize();
      largeRadius_ah /= latDatLBM.GetVoxelSize();
      globalPosition -= latDatLBM.GetOrigin();
      globalPosition /= latDatLBM.GetVoxelSize();

      // updating position with zero velocity and zero body force is necessary
      // because of the side-effect that sets owner rank from the new position
      ownerRank = BIG_NUMBER2;
      velocity *= 0.0;
      bodyForces *= 0.0;
      UpdatePosition(latDatLBM);

      OutputInformation();
      if (log::Logger::ShouldDisplay<log::Debug>())
        log::Logger::Log<log::Debug, log::OnePerCore>(
          "In colloids::Particle::ctor, id: %i, a0: %g, ah: %g, position: {%g,%g,%g}\n",
          particleId, smallRadius_a0, largeRadius_ah,
          globalPosition.x, globalPosition.y, globalPosition.z);
    }

    const bool Particle::operator<(const Particle& other) const
    {
      // ORDER BY isLocal, ownerRank, particleId
      if (ownerRank == other.ownerRank)
        return (particleId < other.particleId);
      else if (ownerRank == topology::NetworkTopology::Instance()->GetLocalRank())
        return true;
      else if (other.ownerRank == topology::NetworkTopology::Instance()->GetLocalRank())
        return false;
      else
        return (ownerRank < other.ownerRank);
    }

    const bool Particle::IsOwnerRankKnown(std::map<proc_t, std::pair<unsigned int, unsigned int> > map) const
    {
      return map.count(ownerRank)>0;
    }

    const void Particle::OutputInformation() const
    {
        log::Logger::Log<log::Info, log::OnePerCore>(
          "In colloids::Particle::OutputInformation, id: %i, owner: %i, position: {%g,%g,%g}, velocity: {%g,%g,%g}, bodyForces: {%g,%g,%g}\n",
          particleId, ownerRank, globalPosition.x, globalPosition.y, globalPosition.z,
          velocity.x, velocity.y, velocity.z, bodyForces.x, bodyForces.y, bodyForces.z);
    }

    const void Particle::UpdatePosition(const geometry::LatticeData& latDatLBM)
    {
      if (log::Logger::ShouldDisplay<log::Debug>())
        log::Logger::Log<log::Debug, log::OnePerCore>(
          "In colloids::Particle::UpdatePosition, id: %i,\nposition: {%g,%g,%g}\nvelocity: {%g,%g,%g}\nbodyForces: {%g,%g,%g}\n",
          particleId, globalPosition.x, globalPosition.y, globalPosition.z,
          velocity.x, velocity.y, velocity.z, bodyForces.x, bodyForces.y, bodyForces.z);

      globalPosition += velocity + bodyForces;

      // round the global position of the particle to the nearest site coordinates
      const util::Vector3D<site_t> siteGlobalPosition(
                                     (site_t)(0.5+globalPosition.x),
                                     (site_t)(0.5+globalPosition.y),
                                     (site_t)(0.5+globalPosition.z));

      // convert the site coordinates into a local site index and find owner rank
      proc_t procId = latDatLBM.GetProcIdFromGlobalCoords(siteGlobalPosition);
      isValid = (procId != BIG_NUMBER2);
      if (isValid && (ownerRank != procId))
      {
        log::Logger::Log<log::Info, log::OnePerCore>(
          "Changing owner of particle %i from %i to %i - %s\n",
          particleId, ownerRank, procId, isValid ? "valid" : "INVALID");
        ownerRank = procId;
      }

      if (log::Logger::ShouldDisplay<log::Debug>())
        log::Logger::Log<log::Debug, log::OnePerCore>(
          "In colloids::Particle::UpdatePosition, id: %i, position is now: {%g,%g,%g}\n",
          particleId, globalPosition.x, globalPosition.y, globalPosition.z);
    }

    /** creates a derived MPI datatype that represents a single particle object
     *  note - this data type uses displacements rather than absolute addresses
     *  refer to Example 4.17 on pp114-117 of the MPI specification version 2.2
     *  when you no longer need this type, remember to call MPI::Datatype::Free
     */
    const MPI::Datatype Particle::CreateMpiDatatypeWithPosition() const
    {
      // MPI::Get_address specifies non-const pointers
      // so need a non-const copy of a particle object
      Particle temp(*this);

      // we want a re-usable MPI data type
      // so we need relative displacements
      MPI::Aint baseAddress = MPI::Get_address(&temp);

      // we have chosen to make each block of fields contain a single field
      // so, the number of field blocks is the same as the number of fields
      int numberOfFieldBlocks = 5;

      // and the length of every field block is one
      int lengthOfEachFieldBlock[] = {1, 1, 1, 1, 1};

      // there is no guarantee that the fields will be contiguous, so
      // the displacement of each field must be determined separately
      MPI::Aint displacementOfEachFieldBlock[numberOfFieldBlocks];
      displacementOfEachFieldBlock[0] = MPI::Get_address(&(temp.particleId)) - baseAddress; 
      displacementOfEachFieldBlock[1] = MPI::Get_address(&(temp.ownerRank)) - baseAddress; 
      displacementOfEachFieldBlock[2] = MPI::Get_address(&(temp.smallRadius_a0)) - baseAddress; 
      displacementOfEachFieldBlock[3] = MPI::Get_address(&(temp.largeRadius_ah)) - baseAddress; 
      displacementOfEachFieldBlock[4] = MPI::Get_address(&(temp.globalPosition)) - baseAddress; 

      // the built-in MPI datatype of each field must match the C++ type
      MPI::Datatype datatypeOfEachFieldBlock[] =
        {MPI::UNSIGNED_LONG, MPI::INT, MPI::DOUBLE, MPI::DOUBLE, MpiDataType<util::Vector3D<double> >()};

      // create a first draft of the MPI datatype for a Particle
      // the lower bound and displacements of fields are correct
      // but the extent may not include the whole derived object
      // specifically, we aren't sending velocity and bodyForces
      MPI::Datatype particleType = MPI::Datatype::Create_struct(
        numberOfFieldBlocks,
        lengthOfEachFieldBlock,
        displacementOfEachFieldBlock,
        datatypeOfEachFieldBlock);

      // obtain the current lower bound for the MPI datatype
      MPI::Aint lowerBound, extent;
      particleType.Get_extent(lowerBound, extent);

      // we can determine the actual extent of a Particle object
      // by concatenating two of them, using a contiguous vector
      // and finding the difference between their base addresses
      std::vector<Particle> tempVectorOfParticle(2, temp);
      extent = MPI::Get_address(&(tempVectorOfParticle[1]))
             - MPI::Get_address(&(tempVectorOfParticle[0]));

      // resize the uncommitted first draft MPI datatype
      // with the current lower bound and the new extent
      particleType = particleType.Create_resized(lowerBound, extent);

      // commit the MPI datatype and return it
      particleType.Commit();
      return particleType;
    }

    /** creates a derived MPI datatype that represents a single particle object
     *  note - this data type uses displacements rather than absolute addresses
     *  refer to Example 4.17 on pp114-117 of the MPI specification version 2.2
     *  when you no longer need this type, remember to call MPI::Datatype::Free
     */
    const MPI::Datatype Particle::CreateMpiDatatypeWithVelocity() const
    {
      // MPI::Get_address specifies non-const pointers
      // so need a non-const copy of a particle object
      Particle temp(*this);

      // we want a re-usable MPI data type
      // so we need relative displacements
      MPI::Aint baseAddress = MPI::Get_address(&temp);

      // we have chosen to make each block of fields contain a single field
      // so, the number of field blocks is the same as the number of fields
      int numberOfFieldBlocks = 2;

      // and the length of every field block is one
      int lengthOfEachFieldBlock[] = {1, 1};

      // there is no guarantee that the fields will be contiguous, so
      // the displacement of each field must be determined separately
      MPI::Aint displacementOfEachFieldBlock[numberOfFieldBlocks];
      displacementOfEachFieldBlock[0] = MPI::Get_address(&(temp.particleId)) - baseAddress; 
      displacementOfEachFieldBlock[1] = MPI::Get_address(&(temp.velocity)) - baseAddress; 

      // the built-in MPI datatype of each field must match the C++ type
      MPI::Datatype datatypeOfEachFieldBlock[] =
        {MPI::UNSIGNED_LONG, MpiDataType<util::Vector3D<double> >()};

      // create a first draft of the MPI datatype for a Particle
      // the lower bound and displacements of fields are correct
      // but the extent may not include the whole derived object
      // specifically, we aren't sending velocity and bodyForces
      MPI::Datatype particleType = MPI::Datatype::Create_struct(
        numberOfFieldBlocks,
        lengthOfEachFieldBlock,
        displacementOfEachFieldBlock,
        datatypeOfEachFieldBlock);

      // obtain the current lower bound for the MPI datatype
      MPI::Aint lowerBound, extent;
      particleType.Get_extent(lowerBound, extent);

      // we can determine the actual extent of a Particle object
      // by concatenating two of them, using a contiguous vector
      // and finding the difference between their base addresses
      std::vector<Particle> tempVectorOfParticle(2, temp);
      extent = MPI::Get_address(&(tempVectorOfParticle[1]))
             - MPI::Get_address(&(tempVectorOfParticle[0]));

      // resize the uncommitted first draft MPI datatype
      // with the current lower bound and the new extent
      particleType = particleType.Create_resized(lowerBound, extent);

      // commit the MPI datatype and return it
      particleType.Commit();
      return particleType;
    }

    const void Particle::CalculateBodyForces()
    {
      /** todo: CalculateBodyForces
       *    - create new BodyForces::GetBodyForcesAtGlobalPosition method
       *    - set bodyForces property to value from the BodyForces object
       *    - will require communicaiton to transmit remote contributions
       */
      bodyForces *= 0.0;
    }

    const void Particle::CalculateFeedbackForces() const
    {
      /** todo: CalculateFeedbackForces
       *    - calculate the feedback force on each neighbour lattice site
       *    - set feedback force values into the macroscopic cache object
       *    - will require communication to transmit remote contributions
       */
    }

    /** modified dirac delta function according to Peskin */
    //double diracOperation(util::Vector3D<double> relativePosition)
    DimensionlessQuantity diracOperation(LatticePosition relativePosition)
    {
      DimensionlessQuantity delta = 1.0;
      for (int xyz=0;xyz<3;xyz++)
      {
        LatticeDistance rmod = fabs(relativePosition[xyz]);

        if (rmod <= 1.0)
          delta *= 0.125*(3.0 - 2.0*rmod + sqrt(1.0 + 4.0*rmod - 4.0*rmod*rmod));
        else if (rmod <= 2.0)
          delta *= 0.125*(5.0 - 2.0*rmod - sqrt(-7.0 + 12.0*rmod  - 4.0*rmod*rmod));
        else
          delta = 0.0;

      }
      return delta;
    }

    const void Particle::InterpolateFluidVelocity(
                           const geometry::LatticeData& latDatLBM,
                           const lb::MacroscopicPropertyCache& propertyCache)
    {
      /** InterpolateFluidVelocity
       *    For each local neighbour lattice site
       *    - get velocity for the site from the macroscopic cache object
       *    - calculate site's contribution to the velocity interpolation
       *    - increment particle velocity property with this contribution
       *    - will require communication to transmit remote contributions
       */

      if (log::Logger::ShouldDisplay<log::Debug>())
        log::Logger::Log<log::Debug, log::OnePerCore>(
          "In colloids::Particle::InterpolateFluidVelocity, id: %i, position: {%g,%g,%g}\n",
          particleId, globalPosition.x, globalPosition.y, globalPosition.z);

      velocity *= 0.0;
      // determine the global coordinates of the next neighbour site:
      // nested loop - x, y, z directions semi-open interval [-2, +2)
      for (site_t x = ((site_t)globalPosition.x)-1;
                  x < ((site_t)globalPosition.x)+3; x++)
        for (site_t y = ((site_t)globalPosition.y)-1;
                    y < ((site_t)globalPosition.y)+3; y++)
          for (site_t z = ((site_t)globalPosition.z)-1;
                      z < ((site_t)globalPosition.z)+3; z++)
          {
            const util::Vector3D<site_t> siteGlobalPosition(x, y, z);

            // convert the global coordinates of the site into a local site index
            proc_t procId;
            site_t siteId;
            bool isSiteValid = latDatLBM.GetContiguousSiteId(siteGlobalPosition, procId, siteId);
            bool isSiteLocal = (procId == topology::NetworkTopology::Instance()->GetLocalRank());

            if (log::Logger::ShouldDisplay<log::Debug>())
            {
              log::Logger::Log<log::Debug, log::OnePerCore>("WAIT A MINUTE 1 ...\n");

              int globalStatus = 0; // valid && local fluid
              int blockStatus = 0; // valid && non-empty
              int siteStatus = 0; // valid && fluid

              util::Vector3D<site_t> blockCoords, localSiteCoords;
              latDatLBM.GetBlockAndLocalSiteCoords(siteGlobalPosition, blockCoords, localSiteCoords);
              log::Logger::Log<log::Debug, log::OnePerCore>("WAIT A MINUTE 2 ...\n");
              if (!latDatLBM.IsValidBlock(blockCoords))
                blockStatus = 1; // invalid - out of range
              else if (latDatLBM.GetBlock(latDatLBM.GetBlockIdFromBlockCoords(blockCoords)).IsEmpty())
                blockStatus = 2; // empty - entire block is solid

              log::Logger::Log<log::Debug, log::OnePerCore>("WAIT A MINUTE 3 ...\n");

              if (!latDatLBM.IsValidLatticeSite(localSiteCoords))
                siteStatus = 1; // invalid - out of range
              else if (blockStatus == 2)
                siteStatus = 2; // deemed solid because block is empty
              else if (blockStatus == 0)
                if (latDatLBM.GetBlock(latDatLBM.GetBlockIdFromBlockCoords(blockCoords)).GetProcessorRankForSite(latDatLBM.GetLocalSiteIdFromLocalSiteCoords(localSiteCoords)) == BIG_NUMBER2)
                  siteStatus = 3; // individual site is not simulated, i.e. must be solid
                else if (latDatLBM.GetBlock(latDatLBM.GetBlockIdFromBlockCoords(blockCoords)).GetProcessorRankForSite(latDatLBM.GetLocalSiteIdFromLocalSiteCoords(localSiteCoords)) == topology::NetworkTopology::Instance()->GetLocalRank())
                  if (latDatLBM.GetBlock(latDatLBM.GetBlockIdFromBlockCoords(blockCoords)).SiteIsSolid(latDatLBM.GetLocalSiteIdFromLocalSiteCoords(localSiteCoords)))
                  siteStatus = 4; // individual site is local but solid (should not happen?)

              log::Logger::Log<log::Debug, log::OnePerCore>("WAIT A MINUTE 4 ...\n");

              if (blockStatus == 1 || (siteStatus == 1))
                globalStatus = 1; // invalid - out of range
              else if (blockStatus == 2 || (siteStatus == 3) | (siteStatus == 4))
                globalStatus = 2; // solid
              else if (latDatLBM.GetBlock(latDatLBM.GetBlockIdFromBlockCoords(blockCoords)).GetProcessorRankForSite(latDatLBM.GetLocalSiteIdFromLocalSiteCoords(localSiteCoords)) != topology::NetworkTopology::Instance()->GetLocalRank())
                globalStatus = 3; // remote fluid

              log::Logger::Log<log::Debug, log::OnePerCore>("ABOUT TO DO STUFF\n");
/*
              log::Logger::Log<log::Debug, log::OnePerCore>(
                "In colloids::Particle::InterpolateFluidVelocity, particleId: %i, siteGlobalPosition: {%i,%i,%i} - %s, blockCoords: {%i,%i,%i} - %s, localSiteCoords: {%i,%i,%i} - %s\n",
                particleId, siteGlobalPosition.x, siteGlobalPosition.y, siteGlobalPosition.z,
                globalStatus == 0 ? "valid && local fluid" :
                (globalStatus == 1 ? "invalid - out-of-range" :
                (globalStatus == 2 ? "solid" :
                (globalStatus == 3 ? "remote fluid" : "unknown status"))),
                particleId, blockCoords.x, blockCoords.y, blockCoords.z,
                blockStatus == 0 ? "valid && non-empty" :
                (blockStatus == 1 ? "invalid - out-of-range" :
                (blockStatus == 2 ? "empty - entire block is solid" : "unknown status")),
                particleId, localSiteCoords.x, localSiteCoords.y, localSiteCoords.z,
                siteStatus == 0 ? "valid && fluid" :
                (siteStatus == 1 ? "invalid - out of range" :
                (siteStatus == 2 ? "deemed solid because block is empty" :
                (siteStatus == 3 ? "individual site is solid" : "unknown status"))));
*/
//*
              log::Logger::Log<log::Debug, log::OnePerCore>(
                "In colloids::Particle::InterpolateFluidVelocity, particleId: %i, siteGlobalPosition: {%i,%i,%i} - %s \n",
                particleId, siteGlobalPosition.x, siteGlobalPosition.y, siteGlobalPosition.z,
                globalStatus == 0 ? "valid && local fluid" :
                (globalStatus == 1 ? "invalid - out-of-range" :
                (globalStatus == 2 ? "solid" :
                (globalStatus == 3 ? "remote fluid" : "unknown status"))));
              log::Logger::Log<log::Debug, log::OnePerCore>("ABOUT TO DO MORE STUFF\n");
              log::Logger::Log<log::Debug, log::OnePerCore>(
                "In colloids::Particle::InterpolateFluidVelocity, particleId: %i, blockCoords: {%i,%i,%i} - %s\n",
                particleId, blockCoords.x, blockCoords.y, blockCoords.z,
                blockStatus == 0 ? "valid && non-empty" :
                (blockStatus == 1 ? "invalid - out-of-range" :
                (blockStatus == 2 ? "empty - entire block is solid" : "unknown status")));
              log::Logger::Log<log::Debug, log::OnePerCore>("ABOUT TO DO EVEN MORE STUFF\n");
              log::Logger::Log<log::Debug, log::OnePerCore>(
                "In colloids::Particle::InterpolateFluidVelocity, particleId: %i, localSiteCoords: {%i,%i,%i} - %s\n",
                particleId, localSiteCoords.x, localSiteCoords.y, localSiteCoords.z,
                siteStatus == 0 ? "valid && fluid" :
                (siteStatus == 1 ? "invalid - out of range" :
                (siteStatus == 2 ? "deemed solid because block is empty" :
                (siteStatus == 3 ? "deemed solid because not simulated" :
                (siteStatus == 4 ? "local but solid (should not happen?)" : "unknown status")))));
//*/
              if (siteStatus == 0)
              log::Logger::Log<log::Debug, log::OnePerCore>("ABOUT TO DO ONE LAST THING\n");
              if (siteStatus == 0)
              log::Logger::Log<log::Debug, log::OnePerCore>(
                "In colloids::Particle::InterpolateFluidVelocity, particleId: %i, procIdForSite: %i\n",
                particleId, latDatLBM.GetBlock(latDatLBM.GetBlockIdFromBlockCoords(blockCoords)).GetProcessorRankForSite(latDatLBM.GetLocalSiteIdFromLocalSiteCoords(localSiteCoords)));
              log::Logger::Log<log::Debug, log::OnePerCore>("DONE LOTS OF STUFF\n");

              if (!isSiteValid || !isSiteLocal)
                if (globalStatus == 0)
                {
                  log::Logger::Log<log::Debug, log::OnePerCore>("MAJOR PROBLEM! isValid: %s, isLocal: %s, procIdForSite: %i, localRank: %i\n", isSiteValid ? "T": "F", isSiteLocal ? "T" : "F",
                    latDatLBM.GetBlock(latDatLBM.GetBlockIdFromBlockCoords(blockCoords)).GetProcessorRankForSite(latDatLBM.GetLocalSiteIdFromLocalSiteCoords(localSiteCoords)),
                    topology::NetworkTopology::Instance()->GetLocalRank());
                }
            }

            /** TODO: implement communication for remote fluid sites */
            /** TODO: implement boundary conditions for invalid/solid sites */
            if (!isSiteValid || !isSiteLocal)
              continue;

            // read value of velocity for site index from macroscopic cache
            // TODO: should be LatticeVelocity == Vector3D<LatticeSpeed> (fix as part of #437)
            util::Vector3D<double> siteFluidVelocity = propertyCache.velocityCache.Get(siteId);

            // calculate term of the interpolation sum
            LatticePosition relativePosition(siteGlobalPosition);
            relativePosition -= globalPosition;
            // TODO: should be LatticeVelocity == Vector3D<LatticeSpeed> (fix as part of #437)
            util::Vector3D<double> partialInterpolation = siteFluidVelocity *
              diracOperation(relativePosition);

            // accumulate each term of the interpolation
            velocity += partialInterpolation;

            if (log::Logger::ShouldDisplay<log::Debug>())
              log::Logger::Log<log::Debug, log::OnePerCore>(
                "In colloids::Particle::InterpolateFluidVelocity, particleId: %i, siteIndex: %i, fluidVelocity: {%g,%g,%g}, partialInterpolation: {%g,%g,%g}, velocitySoFar: {%g,%g,%g}\n",
                particleId, siteId, siteFluidVelocity.x, siteFluidVelocity.y, siteFluidVelocity.z,
                partialInterpolation.x, partialInterpolation.y, partialInterpolation.z,
                velocity.x, velocity.y, velocity.z);
          }
    }

  }
}
