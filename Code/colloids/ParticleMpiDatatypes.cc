// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#include "colloids/Particle.h"
//#include "units.h"

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
      int numberOfFieldBlocks = 6;

      // and the length of every field block is one
      int lengthOfEachFieldBlock[] = {1, 1, 1, 1, 1, 1};

      // there is no guarantee that the fields will be contiguous, so
      // the displacement of each field must be determined separately
      MPI::Aint displacementOfEachFieldBlock[numberOfFieldBlocks];
      displacementOfEachFieldBlock[0] = MPI::Get_address(&(temp.particleId)) - baseAddress; 
      displacementOfEachFieldBlock[1] = MPI::Get_address(&(temp.ownerRank)) - baseAddress; 
      displacementOfEachFieldBlock[2] = MPI::Get_address(&(temp.smallRadius_a0)) - baseAddress; 
      displacementOfEachFieldBlock[3] = MPI::Get_address(&(temp.largeRadius_ah)) - baseAddress; 
      displacementOfEachFieldBlock[4] = MPI::Get_address(&(temp.mass)) - baseAddress; 
      displacementOfEachFieldBlock[5] = MPI::Get_address(&(temp.globalPosition)) - baseAddress; 

      // the built-in MPI datatype of each field must match the C++ type
      MPI::Datatype datatypeOfEachFieldBlock[] =
        {MPI::UNSIGNED_LONG, MPI::INT, MPI::DOUBLE, MPI::DOUBLE,
         MPI::DOUBLE, MpiDataType<util::Vector3D<double> >()};

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

  }
}
