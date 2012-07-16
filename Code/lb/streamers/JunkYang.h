// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_LB_STREAMERS_JUNKYANG_H
#define HEMELB_LB_STREAMERS_JUNKYANG_H

#include "lb/kernels/BaseKernel.h"
#include "lb/streamers/BaseStreamer.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace hemelb
{
  namespace lb
  {
    namespace streamers
    {
      using namespace boost::numeric;

      /**
       * This class implements the Junk&Yang boundary condition as described in
       *
       * M. Junk and Z. Yang "One-point boundary condition for the lattice Boltzmann method", Phys Rev E 72 (2005)
       */
      template<typename CollisionImpl>
      class JunkYang : public BaseStreamer<JunkYang<CollisionImpl> >
      {
        public:

          typedef CollisionImpl CollisionType;

          //! @todo: no default constructor for ublas::permutation_matrix. Rewrite it as a vector of pointers and do the allocation/deallocation
          JunkYang(kernels::InitParams& initParams) :
              collider(initParams), THETA(1.0), kernelFirstSite(initParams.firstSite), kernelSiteCount(initParams.siteCount), latticeData(*initParams.latDat), incomingVelocities(kernelSiteCount), outgoingVelocities(kernelSiteCount), kMatrices(kernelSiteCount), lMatrices(kernelSiteCount), luPermutationMatrices(kernelSiteCount,
                                                                                                                                                                                                                                                                                                              ublas::permutation_matrix<
                                                                                                                                                                                                                                                                                                                  size_t>(LatticeType::NUMVECTORS)), fPostCollision(kernelSiteCount), fPostCollisionInverseDir(kernelSiteCount), fOld(kernelSiteCount)
          {
            AssembleMatrices();
            FactoriseLMatrices();
          }

          template<bool tDoRayTracing>
          inline void DoStreamAndCollide(const site_t firstIndex,
                                         const site_t siteCount,
                                         const LbmParameters* lbmParams,
                                         geometry::LatticeData* latticeData,
                                         lb::MacroscopicPropertyCache& propertyCache)
          {
            assert(firstIndex == kernelFirstSite);
            assert(siteCount == kernelSiteCount);

            for (site_t siteIndex = firstIndex; siteIndex < (firstIndex + siteCount); siteIndex++)
            {
              geometry::Site site = latticeData->GetSite(siteIndex);

              distribn_t *distribution = site.GetFOld<LatticeType>();

              kernels::HydroVars<typename CollisionType::CKernel> hydroVars(distribution);

              ///< @todo #126 This value of tau will be updated by some kernels within the collider code (e.g. LBGKNN). It would be nicer if tau is handled in a single place.
              hydroVars.tau = lbmParams->GetTau();

              // In the first step, we stream and collide as we would for the SimpleCollideAndStream
              // streamer.
              collider.CalculatePreCollision(hydroVars, site);
              collider.Collide(lbmParams, hydroVars);

              for (unsigned int direction = 0; direction < LatticeType::NUMVECTORS; direction++)
              {
                if (!site.HasBoundary(direction))
                {
                  distribution[direction] = hydroVars.GetFPostCollision()[direction];

                  * (latticeData->GetFNew(site.GetStreamedIndex<LatticeType>(direction))) =
                      hydroVars.GetFPostCollision()[direction];
                }
              }

              // The following for loops prepare the data required by DoPostStep
              unsigned siteLocalIndex = siteIndex - firstIndex;
              fPostCollisionInverseDir[siteLocalIndex].resize(incomingVelocities[siteLocalIndex].size());

              unsigned index = 0;
              for (std::set<Direction>::const_iterator incomingVelocityIter =
                  incomingVelocities[siteLocalIndex].begin();
                  incomingVelocityIter != incomingVelocities[siteLocalIndex].end(); ++incomingVelocityIter, ++index)
              {
                int inverseDirection = LatticeType::INVERSEDIRECTIONS[*incomingVelocityIter];

                fPostCollision[siteLocalIndex](index) = hydroVars.GetFPostCollision()[*incomingVelocityIter];
                fPostCollisionInverseDir[siteLocalIndex](index) = hydroVars.GetFPostCollision()[inverseDirection];
                fOld[siteLocalIndex](index) = distribution[*incomingVelocityIter];
              }

              for (std::set<Direction>::const_iterator outgoingVelocityIter =
                  outgoingVelocities[siteLocalIndex].begin();
                  outgoingVelocityIter != outgoingVelocities[siteLocalIndex].end(); ++outgoingVelocityIter, ++index)
              {
                fPostCollision[siteLocalIndex](index) = hydroVars.GetFPostCollision()[*outgoingVelocityIter];
                fOld[siteLocalIndex](index) = distribution[*outgoingVelocityIter];
              }

              BaseStreamer<JunkYang>::template UpdateMinsAndMaxes<tDoRayTracing>(hydroVars.v_x,
                                                                                 hydroVars.v_y,
                                                                                 hydroVars.v_z,
                                                                                 site,
                                                                                 hydroVars.GetFNeq().f,
                                                                                 hydroVars.density,
                                                                                 hydroVars.tau,
                                                                                 lbmParams,
                                                                                 propertyCache);
            }
          }

          template<bool tDoRayTracing>
          inline void DoPostStep(const site_t firstIndex,
                                 const site_t siteCount,
                                 const LbmParameters* lbmParams,
                                 geometry::LatticeData* latticeData,
                                 lb::MacroscopicPropertyCache& propertyCache)
          {
            assert(firstIndex == kernelFirstSite);
            assert(siteCount == kernelSiteCount);

            for (site_t siteIndex = firstIndex; siteIndex < (firstIndex + siteCount); siteIndex++)
            {
              unsigned siteLocalIndex = siteIndex - firstIndex;

              ublas::vector < distribn_t > rVector;
              AssembleRVector(siteLocalIndex, rVector);

              // assemble RHS
              ublas::vector < distribn_t > systemRHS = fPostCollisionInverseDir[siteLocalIndex] - rVector;

              // lu_substitue will overwrite the RHS with the solution
              ublas::vector < distribn_t > &systemSolution = systemRHS;
              lu_substitute(lMatrices[siteLocalIndex], luPermutationMatrices[siteLocalIndex], systemSolution);

              // Update the distribution function for incoming velocities with the solution of the linear system
              unsigned index = 0;
              for (std::set<Direction>::const_iterator incomingVelocityIter =
                  incomingVelocities[siteLocalIndex].begin();
                  incomingVelocityIter != incomingVelocities[siteLocalIndex].end(); ++incomingVelocityIter, ++index)
              {
                * (latticeData->GetFNew(siteIndex * LatticeType::NUMVECTORS + *incomingVelocityIter)) =
                    systemSolution[index];
              }
            }
          }

          inline void DoReset(kernels::InitParams* init)
          {
            collider.Reset(init);
          }

        private:

          typedef typename CollisionType::CKernel::LatticeType LatticeType;
          CollisionType collider;

          //! Coordinate index arbitrarily chosen in the paper
          static const unsigned ALPHA = 3u;
          //! theta constant in the theta-method used for interpolation (0 for fully explicit, 1 for fully implicit)
          const distribn_t THETA;

          //! Index of the first local site assigned to this kernel
          site_t kernelFirstSite;
          //! Number of local sites assigned to this kernel
          site_t kernelSiteCount;
          //! Reference to the lattice object used for initialisation
          const geometry::LatticeData& latticeData;

          //! Set of incoming velocities (those with an inverse direction crossing a wall boundary)
          std::vector<std::set<Direction> > incomingVelocities;
          //! Set of outgoing velocities (the complement of incomingVelocities)
          std::vector<std::set<Direction> > outgoingVelocities;

          //! kernelSiteCount long vector containing the K matrices used to assemble both the left-hand- and the right-hand-side of the linear systems for each site
          std::vector<ublas::matrix<distribn_t> > kMatrices;
          //! kernelSiteCount long vector containing the linear system left-hand-sides for each site
          std::vector<ublas::matrix<distribn_t> > lMatrices;
          //! kernelSiteCount long vector containing the permutations generated by the LU factorisation of each linear system
          std::vector<ublas::permutation_matrix<std::size_t> > luPermutationMatrices;

          //! Postcollision distribution with incoming/outgoing ordering.
          std::vector<ublas::c_vector<distribn_t, LatticeType::NUMVECTORS> > fPostCollision;
          //! Postcollision distribution for the inverse of the incoming set of velocities.
          std::vector<ublas::vector<distribn_t> > fPostCollisionInverseDir;
          //! Particle distribution in the previous time step with incoming/outgoing ordering.
          std::vector<ublas::c_vector<distribn_t, LatticeType::NUMVECTORS> > fOld;

          /**
           * Assemble all the required matrices for all the sites assigned to this kernel
           */
          inline void AssembleMatrices()
          {
            for (site_t siteLocalIndex = 0; siteLocalIndex < kernelSiteCount; siteLocalIndex++)
            {
              ConstructVelocitySets(siteLocalIndex);
              AssembleKMatrix(siteLocalIndex);
              AssembleLMatrix(siteLocalIndex);
            }
          }

          /**
           * Construct the incoming/outgoing velocity sets for site siteLocalIndex
           * @param siteLocalIndex site local index (indexed from 0 for the current kernel)
           */
          inline void ConstructVelocitySets(site_t siteLocalIndex)
          {
            geometry::ConstSite site = latticeData.GetSite(kernelFirstSite + siteLocalIndex);

            for (Direction direction = 0; direction < LatticeType::NUMVECTORS; direction++)
            {
              int inverseDirection = LatticeType::INVERSEDIRECTIONS[direction];
              if (site.HasBoundary(inverseDirection))
              {
                incomingVelocities[siteLocalIndex].insert(direction);
              }
              else
              {
                outgoingVelocities[siteLocalIndex].insert(direction);
              }
            }
          }

          /**
           * Assemble the K matrix for site siteLocalIndex.
           *
           * K is a rectangular matrix (num_incoming_vels x LatticeType::NUMVECTORS). Our
           * implementation place the columns corresponding to the set of incoming velocities
           * first followed by outgoing velocities. This makes easier handling the subblocks
           * with uBLAS.
           *
           * @param siteLocalIndex site local index (indexed from 0 for the current kernel)
           */
          inline void AssembleKMatrix(unsigned siteLocalIndex)
          {
            geometry::ConstSite site = latticeData.GetSite(kernelFirstSite + siteLocalIndex);

            unsigned incomingVelsSetSize = incomingVelocities[siteLocalIndex].size();

            kMatrices[siteLocalIndex].resize(incomingVelsSetSize, LatticeType::NUMVECTORS);

            /*
             *  Assemble K(:, 0:num_incoming_velocities-1)
             */
            unsigned rowIndex = 0;
            for (std::set<Direction>::const_iterator incomingVelocityIter = incomingVelocities[siteLocalIndex].begin();
                incomingVelocityIter != incomingVelocities[siteLocalIndex].end(); ++incomingVelocityIter, ++rowIndex)
            {
              // |c_i|^2, where c_i is the i-th velocity vector
              const int rowRowdirectionsInnProd = LatticeType::CX[*incomingVelocityIter]
                  * LatticeType::CX[*incomingVelocityIter]
                  + LatticeType::CY[*incomingVelocityIter] * LatticeType::CY[*incomingVelocityIter]
                  + LatticeType::CZ[*incomingVelocityIter] * LatticeType::CZ[*incomingVelocityIter];

              unsigned columnIndex = 0;
              for (std::set<Direction>::const_iterator incomingVelocityIter2 =
                  incomingVelocities[siteLocalIndex].begin();
                  incomingVelocityIter2 != incomingVelocities[siteLocalIndex].end();
                  ++incomingVelocityIter2, ++columnIndex)
              {
                // |c_i|^2, where c_i is the i-th velocity vector
                const int colColdirectionsInnProd = LatticeType::CX[*incomingVelocityIter2]
                    * LatticeType::CX[*incomingVelocityIter2]
                    + LatticeType::CY[*incomingVelocityIter2] * LatticeType::CY[*incomingVelocityIter2]
                    + LatticeType::CZ[*incomingVelocityIter2] * LatticeType::CZ[*incomingVelocityIter2];

                // (c_i \dot c_j)^2, where c_{i,j} are the {i,j}-th velocity vectors
                const int rowColdirectionsInnProd = abs(LatticeType::CX[*incomingVelocityIter]
                    * LatticeType::CX[*incomingVelocityIter2]
                    + LatticeType::CY[*incomingVelocityIter] * LatticeType::CY[*incomingVelocityIter2]
                    + LatticeType::CZ[*incomingVelocityIter] * LatticeType::CZ[*incomingVelocityIter2]);

                assert(site.template GetWallDistance<LatticeType>(LatticeType::INVERSEDIRECTIONS[(Direction) *incomingVelocityIter])>=0);
                assert(site.template GetWallDistance<LatticeType>(LatticeType::INVERSEDIRECTIONS[(Direction) *incomingVelocityIter])<1);

                kMatrices[siteLocalIndex](rowIndex, columnIndex) =
                    -3.0 / 2.0
                        * (3.0
                            - 6
                                * site.template GetWallDistance<LatticeType>(LatticeType::INVERSEDIRECTIONS[(Direction) *incomingVelocityIter]))
                        * LatticeType::EQMWEIGHTS[*incomingVelocityIter]
                        * ( (rowColdirectionsInnProd * rowColdirectionsInnProd) - (rowRowdirectionsInnProd / 3.0)
                            - LatticeType::discreteVelocityVectors[ALPHA - 1][*incomingVelocityIter]
                                * (colColdirectionsInnProd - (ALPHA / 3.0)));

                assert(fabs(kMatrices[siteLocalIndex](rowIndex, columnIndex)) < 1e3);
              }

              /*
               * Assemble K(:, num_incoming_velocities:num_incoming_velocities+num_outgoing_velocities-1)
               */
              for (std::set<Direction>::const_iterator outgoingVelocityIter =
                  outgoingVelocities[siteLocalIndex].begin();
                  outgoingVelocityIter != outgoingVelocities[siteLocalIndex].end();
                  ++outgoingVelocityIter, ++columnIndex)
              {
                // |c_i|^2, where c_i is the i-th velocity vector
                const int colColdirectionsInnProd = LatticeType::CX[*outgoingVelocityIter]
                    * LatticeType::CX[*outgoingVelocityIter]
                    + LatticeType::CY[*outgoingVelocityIter] * LatticeType::CY[*outgoingVelocityIter]
                    + LatticeType::CZ[*outgoingVelocityIter] * LatticeType::CZ[*outgoingVelocityIter];

                // (c_i \dot c_j)^2, where c_{i,j} are the {i,j}-th velocity vectors
                const int rowColdirectionsInnProd = abs(LatticeType::CX[*incomingVelocityIter]
                    * LatticeType::CX[*outgoingVelocityIter]
                    + LatticeType::CY[*incomingVelocityIter] * LatticeType::CY[*outgoingVelocityIter]
                    + LatticeType::CZ[*incomingVelocityIter] * LatticeType::CZ[*outgoingVelocityIter]);

                assert(site.template GetWallDistance<LatticeType>(LatticeType::INVERSEDIRECTIONS[(Direction) *incomingVelocityIter])>=0);
                assert(site.template GetWallDistance<LatticeType>(LatticeType::INVERSEDIRECTIONS[(Direction) *incomingVelocityIter])<1);

                kMatrices[siteLocalIndex](rowIndex, columnIndex) =
                    -3.0 / 2.0
                        * (3.0
                            - 6
                                * site.template GetWallDistance<LatticeType>(LatticeType::INVERSEDIRECTIONS[(Direction) *incomingVelocityIter]))
                        * LatticeType::EQMWEIGHTS[*incomingVelocityIter]
                        * ( (rowColdirectionsInnProd * rowColdirectionsInnProd) - (rowRowdirectionsInnProd / 3.0)
                            - LatticeType::discreteVelocityVectors[ALPHA - 1][*incomingVelocityIter]
                                * (colColdirectionsInnProd - (ALPHA / 3.0)));

                assert(fabs(kMatrices[siteLocalIndex](rowIndex, columnIndex)) < 1e3);
              }

            }
          }

          /**
           * Assemble the L matrix for for site siteLocalIndex. L is a square matrix (num_incoming_vels x num_incoming_vels)
           *
           * @param siteLocalIndex site local index (indexed from 0 for the current kernel)
           */
          inline void AssembleLMatrix(site_t siteLocalIndex)
          {
            unsigned incomingVelsSetSize = incomingVelocities[siteLocalIndex].size();

            lMatrices[siteLocalIndex] = ublas::identity_matrix<distribn_t>(incomingVelsSetSize)
                + THETA * ublas::subrange(kMatrices[siteLocalIndex], 0, incomingVelsSetSize, 0, incomingVelsSetSize);
          }

          /**
           * Compute the LU factorisation of all the L matrices
           */
          inline void FactoriseLMatrices()
          {
            for (site_t siteLocalIndex = 0; siteLocalIndex < kernelSiteCount; siteLocalIndex++)
            {
              luPermutationMatrices[siteLocalIndex].resize(incomingVelocities[siteLocalIndex].size());
              int ret = lu_factorize(lMatrices[siteLocalIndex], luPermutationMatrices[siteLocalIndex]);
              // If this assertion trips, lMatrices[siteLocalIndex] is singular.
              assert(ret == 0);
            }
          }

          /**
           * Assemble the sigma vector required to assemble the system RHS. We are not including
           * the forcing term used in the paper to drive the flow. This might become necessary
           * for biocolloids.
           *
           * @param siteLocalIndex site local index (indexed from 0 for the current kernel)
           * @param sigmaVector sigma vector
           */
          inline void AssembleSigmaVector(const site_t siteLocalIndex,
                                          ublas::c_vector<distribn_t, LatticeType::NUMVECTORS>& sigmaVector) const
          {
            sigmaVector = fPostCollision[siteLocalIndex] - (1 - THETA) * fOld[siteLocalIndex];
          }

          /**
           * Assemble the r vector required to assemble the system RHS
           *
           * @param siteLocalIndex site local index (indexed from 0 for the current kernel)
           * @param rVector r vector
           */
          inline void AssembleRVector(const site_t siteLocalIndex,
                                      ublas::vector<distribn_t>& rVector) const
          {
            ublas::c_vector < distribn_t, LatticeType::NUMVECTORS > sigmaVector;
            AssembleSigmaVector(siteLocalIndex, sigmaVector);

            unsigned incomingVelsSetSize = incomingVelocities[siteLocalIndex].size();
            unsigned outgoingVelsSetSize = outgoingVelocities[siteLocalIndex].size();

            ublas::vector < distribn_t > fNew(outgoingVelsSetSize);
            site_t siteIndex = kernelFirstSite + siteLocalIndex;

            /*
             *  Assemble a vector with the updated values of the distribution function for the outgoing velocities,
             *  which have already been streamed
             */
            assert(fNew.size() == outgoingVelocities[siteLocalIndex].size());
            unsigned index = 0;
            for (std::set<Direction>::const_iterator outgoingDirIter = outgoingVelocities[siteLocalIndex].begin();
                outgoingDirIter != outgoingVelocities[siteLocalIndex].end(); ++outgoingDirIter, ++index)
            {
              fNew[index] = *latticeData.GetFNew(siteIndex * LatticeType::NUMVECTORS + *outgoingDirIter);
            }

            rVector = THETA
                * ublas::prod(ublas::subrange(kMatrices[siteLocalIndex],
                                              0,
                                              incomingVelsSetSize,
                                              incomingVelsSetSize,
                                              LatticeType::NUMVECTORS),
                              fNew) + ublas::prod(kMatrices[siteLocalIndex], sigmaVector);
          }
      };
    }
  }
}

#endif /* HEMELB_LB_STREAMERS_JUNKYANG_H */
