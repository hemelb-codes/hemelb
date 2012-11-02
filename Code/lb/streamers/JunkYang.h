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

          JunkYang(kernels::InitParams& initParams) :
            collider(initParams), THETA(0.7), latticeData(*initParams.latDat)
          {
            AssembleMatrices();
            FactoriseLMatrices();
          }

          ~JunkYang()
          {
            // Free dynamically allocated permutation matrices.
            for (std::map<site_t, ublas::permutation_matrix<std::size_t>*>::iterator iter =
                luPermutationMatrices.begin(); iter != luPermutationMatrices.end(); ++iter)
            {
              delete iter->second;
            }
          }

          template<bool tDoRayTracing>
          inline void DoStreamAndCollide(const site_t firstIndex,
                                         const site_t siteCount,
                                         const LbmParameters* lbmParams,
                                         geometry::LatticeData* latticeData,
                                         lb::MacroscopicPropertyCache& propertyCache)
          {
            for (site_t siteIndex = firstIndex; siteIndex < (firstIndex + siteCount); siteIndex++)
            {
              assert(latticeData->GetSite(siteIndex).IsEdge());
              assert(lMatrices.find(siteIndex) != lMatrices.end());

              geometry::Site<geometry::LatticeData> site = latticeData->GetSite(siteIndex);

              kernels::HydroVars<typename CollisionType::CKernel> hydroVars(site.GetFOld<LatticeType> ());

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
                  * (latticeData->GetFNew(site.GetStreamedIndex<LatticeType> (direction)))
                      = hydroVars.GetFPostCollision()[direction];
                }
              }

              // The following for loops prepare the data required by DoPostStep
              fPostCollisionInverseDir[siteIndex].resize(incomingVelocities[siteIndex].size());

              unsigned index = 0;
              for (std::set<Direction>::const_iterator incomingVelocityIter = incomingVelocities[siteIndex].begin(); incomingVelocityIter
                  != incomingVelocities[siteIndex].end(); ++incomingVelocityIter, ++index)
              {
                int inverseDirection = LatticeType::INVERSEDIRECTIONS[*incomingVelocityIter];

                fPostCollision[siteIndex](index) = hydroVars.GetFPostCollision()[*incomingVelocityIter];
                fPostCollisionInverseDir[siteIndex](index) = hydroVars.GetFPostCollision()[inverseDirection];
                fOld[siteIndex](index) = site.GetFOld<LatticeType> ()[*incomingVelocityIter];
              }

              for (std::set<Direction>::const_iterator outgoingVelocityIter = outgoingVelocities[siteIndex].begin(); outgoingVelocityIter
                  != outgoingVelocities[siteIndex].end(); ++outgoingVelocityIter, ++index)
              {
                fPostCollision[siteIndex](index) = hydroVars.GetFPostCollision()[*outgoingVelocityIter];
                fOld[siteIndex](index) = site.GetFOld<LatticeType> ()[*outgoingVelocityIter];
              }

              BaseStreamer<JunkYang>::template UpdateMinsAndMaxes<tDoRayTracing>(site,
                                                                                 hydroVars,
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
            for (site_t siteIndex = firstIndex; siteIndex < (firstIndex + siteCount); siteIndex++)
            {
              assert(latticeData->GetSite(siteIndex).IsEdge());
              assert(lMatrices.find(siteIndex) != lMatrices.end());
              assert(luPermutationMatrices.find(siteIndex) != luPermutationMatrices.end());

              ublas::vector < distribn_t > rVector;
              AssembleRVector(siteIndex, rVector);

              // assemble RHS
              ublas::vector < distribn_t > systemRHS = fPostCollisionInverseDir[siteIndex] - rVector;

              // lu_substitue will overwrite the RHS with the solution
              ublas::vector < distribn_t > &systemSolution = systemRHS;
              lu_substitute(lMatrices[siteIndex], *luPermutationMatrices[siteIndex], systemSolution);

              // Update the distribution function for incoming velocities with the solution of the linear system
              unsigned index = 0;
              for (std::set<Direction>::const_iterator incomingVelocityIter = incomingVelocities[siteIndex].begin(); incomingVelocityIter
                  != incomingVelocities[siteIndex].end(); ++incomingVelocityIter, ++index)
              {
                * (latticeData->GetFNew(siteIndex * LatticeType::NUMVECTORS + *incomingVelocityIter))
                    = systemSolution[index];
              }
            }
          }

        private:

          typedef typename CollisionType::CKernel::LatticeType LatticeType;
          CollisionType collider;

          //! Problem dimension (2D, 3D)
          static const unsigned DIMENSION = 3U;
          //! Vector coordinate arbitrarily chosen in the paper
          static const unsigned ALPHA = DIMENSION - 1;
          //! theta constant in the theta-method used for interpolation (0 for fully explicit, 1 for fully implicit)
          const distribn_t THETA;

          //! Reference to the lattice object used for initialisation
          const geometry::LatticeData& latticeData;

          //! Sets of incoming velocities (those with an inverse direction crossing a wall boundary)
          std::map<site_t, std::set<Direction> > incomingVelocities;
          //! Sets of outgoing velocities (the complement of incomingVelocities)
          std::map<site_t, std::set<Direction> > outgoingVelocities;

          //! Map containing the K matrices used to assemble both the left-hand- and the right-hand-side of the linear systems for each site
          std::map<site_t, ublas::matrix<distribn_t> > kMatrices;
          //! Map containing the linear system left-hand-sides for each site
          std::map<site_t, ublas::matrix<distribn_t> > lMatrices;
          /**
           * Map containing the permutations generated by the LU factorisation of each linear system. The permutaation matrices are dinamically allocated
           * because there isn't a default constructor for ublas::permutation_matrix (unlike other ublas data structures) and the size is not known until
           * run time.
           */
          std::map<site_t, ublas::permutation_matrix<std::size_t>*> luPermutationMatrices;

          //! Map of postcollision distributions with incoming/outgoing ordering.
          std::map<site_t, ublas::c_vector<distribn_t, LatticeType::NUMVECTORS> > fPostCollision;
          //! Map of postcollision distributions for the inverse directions of those in the incoming set of velocities.
          std::map<site_t, ublas::vector<distribn_t> > fPostCollisionInverseDir;
          //! Map of distributions in the previous time step with incoming/outgoing ordering.
          std::map<site_t, ublas::c_vector<distribn_t, LatticeType::NUMVECTORS> > fOld;

          /**
           * Assemble all the required matrices for all the sites assigned to this kernel
           */
          inline void AssembleMatrices()
          {
            for (site_t contiguousSiteIndex = 0; contiguousSiteIndex < latticeData.GetLocalFluidSiteCount(); ++contiguousSiteIndex)
            {
              geometry::Site<const geometry::LatticeData> localSite = latticeData.GetSite(contiguousSiteIndex);

              // Ignore ones that aren't edges;
              //! @todo: We should also be ignoring iolets
              if (!localSite.IsEdge())
              {
                continue;
              }

              ConstructVelocitySets(contiguousSiteIndex);
              AssembleKMatrix(contiguousSiteIndex);
              AssembleLMatrix(contiguousSiteIndex);
            }
          }

          /**
           * Construct the incoming/outgoing velocity sets for site siteLocalIndex
           *
           * @param contiguousSiteIndex Contiguous site index (for this core)
           */
          inline void ConstructVelocitySets(site_t contiguousSiteIndex)
          {
            geometry::Site<const geometry::LatticeData> site = latticeData.GetSite(contiguousSiteIndex);

            for (Direction direction = 0; direction < LatticeType::NUMVECTORS; direction++)
            {
              int inverseDirection = LatticeType::INVERSEDIRECTIONS[direction];
              if (site.HasBoundary(inverseDirection))
              {
                incomingVelocities[contiguousSiteIndex].insert(direction);
              }
              else
              {
                outgoingVelocities[contiguousSiteIndex].insert(direction);
              }
            }
          }

          /**
           * Assemble the K matrix for site contiguousSiteIndex.
           *
           * K is a rectangular matrix (num_incoming_vels x LatticeType::NUMVECTORS). Our
           * implementation places the columns corresponding to the set of incoming velocities
           * first followed by those corresponding to outgoing velocities. This makes easier
           * handling the subblocks with uBLAS.
           *
           * @param contiguousSiteIndex Contiguous site index (for this core)
           */
          inline void AssembleKMatrix(site_t contiguousSiteIndex)
          {
            geometry::Site<const geometry::LatticeData> site = latticeData.GetSite(contiguousSiteIndex);

            unsigned incomingVelsSetSize = incomingVelocities[contiguousSiteIndex].size();

            kMatrices[contiguousSiteIndex].resize(incomingVelsSetSize, LatticeType::NUMVECTORS);

            /*
             *  Assemble K(:, 0:num_incoming_velocities-1)
             */
            unsigned rowIndex = 0;
            for (std::set<Direction>::const_iterator rowIndexIncomingVelocity =
                incomingVelocities[contiguousSiteIndex].begin(); rowIndexIncomingVelocity
                != incomingVelocities[contiguousSiteIndex].end(); ++rowIndexIncomingVelocity, ++rowIndex)
            {
              // |c_i|^2, where c_i is the i-th velocity vector
              const int rowRowdirectionsInnProd = LatticeType::CX[*rowIndexIncomingVelocity]
                  * LatticeType::CX[*rowIndexIncomingVelocity] + LatticeType::CY[*rowIndexIncomingVelocity]
                  * LatticeType::CY[*rowIndexIncomingVelocity] + LatticeType::CZ[*rowIndexIncomingVelocity]
                  * LatticeType::CZ[*rowIndexIncomingVelocity];

              unsigned columnIndex = 0;
              for (std::set<Direction>::const_iterator columnIndexIncomingVelocity =
                  incomingVelocities[contiguousSiteIndex].begin(); columnIndexIncomingVelocity
                  != incomingVelocities[contiguousSiteIndex].end(); ++columnIndexIncomingVelocity, ++columnIndex)
              {
                // |c_i|^2, where c_i is the i-th velocity vector
                const int colColdirectionsInnProd = LatticeType::CX[*columnIndexIncomingVelocity]
                    * LatticeType::CX[*columnIndexIncomingVelocity] + LatticeType::CY[*columnIndexIncomingVelocity]
                    * LatticeType::CY[*columnIndexIncomingVelocity] + LatticeType::CZ[*columnIndexIncomingVelocity]
                    * LatticeType::CZ[*columnIndexIncomingVelocity];

                // (c_i \dot c_j)^2, where c_{i,j} are the {i,j}-th velocity vectors
                const int rowColdirectionsInnProd = LatticeType::CX[*rowIndexIncomingVelocity]
                    * LatticeType::CX[*columnIndexIncomingVelocity] + LatticeType::CY[*rowIndexIncomingVelocity]
                    * LatticeType::CY[*columnIndexIncomingVelocity] + LatticeType::CZ[*rowIndexIncomingVelocity]
                    * LatticeType::CZ[*columnIndexIncomingVelocity];

                assert(site.template GetWallDistance<LatticeType> (LatticeType::INVERSEDIRECTIONS[(Direction) *rowIndexIncomingVelocity])
                    >= 0);
                assert(site.template GetWallDistance<LatticeType> (LatticeType::INVERSEDIRECTIONS[(Direction) *rowIndexIncomingVelocity])
                    < 1);

                kMatrices[contiguousSiteIndex](rowIndex, columnIndex)
                    = (-3.0 / 2.0)
                        * (3.0
                            - 6
                                * site.template GetWallDistance<LatticeType> (LatticeType::INVERSEDIRECTIONS[(Direction) *rowIndexIncomingVelocity]))
                        * LatticeType::EQMWEIGHTS[*rowIndexIncomingVelocity] * ( (rowColdirectionsInnProd
                        * rowColdirectionsInnProd) - (rowRowdirectionsInnProd / 3.0)
                        - LatticeType::discreteVelocityVectors[ALPHA][*rowIndexIncomingVelocity]
                            * LatticeType::discreteVelocityVectors[ALPHA][*rowIndexIncomingVelocity]
                            * (colColdirectionsInnProd - (DIMENSION / 3.0)));

                assert(fabs(kMatrices[contiguousSiteIndex](rowIndex, columnIndex)) < 1e3);
              }

              /*
               * Assemble K(:, num_incoming_velocities:num_incoming_velocities+num_outgoing_velocities-1)
               */
              for (std::set<Direction>::const_iterator columnIndexOutgoingVelocity =
                  outgoingVelocities[contiguousSiteIndex].begin(); columnIndexOutgoingVelocity
                  != outgoingVelocities[contiguousSiteIndex].end(); ++columnIndexOutgoingVelocity, ++columnIndex)
              {
                // |c_i|^2, where c_i is the i-th velocity vector
                const int colColdirectionsInnProd = LatticeType::CX[*columnIndexOutgoingVelocity]
                    * LatticeType::CX[*columnIndexOutgoingVelocity] + LatticeType::CY[*columnIndexOutgoingVelocity]
                    * LatticeType::CY[*columnIndexOutgoingVelocity] + LatticeType::CZ[*columnIndexOutgoingVelocity]
                    * LatticeType::CZ[*columnIndexOutgoingVelocity];

                // (c_i \dot c_j)^2, where c_{i,j} are the {i,j}-th velocity vectors
                const int rowColdirectionsInnProd = LatticeType::CX[*rowIndexIncomingVelocity]
                    * LatticeType::CX[*columnIndexOutgoingVelocity] + LatticeType::CY[*rowIndexIncomingVelocity]
                    * LatticeType::CY[*columnIndexOutgoingVelocity] + LatticeType::CZ[*rowIndexIncomingVelocity]
                    * LatticeType::CZ[*columnIndexOutgoingVelocity];

                assert(site.template GetWallDistance<LatticeType> (LatticeType::INVERSEDIRECTIONS[(Direction) *rowIndexIncomingVelocity])
                    >= 0);
                assert(site.template GetWallDistance<LatticeType> (LatticeType::INVERSEDIRECTIONS[(Direction) *rowIndexIncomingVelocity])
                    < 1);

                kMatrices[contiguousSiteIndex](rowIndex, columnIndex)
                    = (-3.0 / 2.0)
                        * (3.0
                            - 6
                                * site.template GetWallDistance<LatticeType> (LatticeType::INVERSEDIRECTIONS[(Direction) *rowIndexIncomingVelocity]))
                        * LatticeType::EQMWEIGHTS[*rowIndexIncomingVelocity] * ( (rowColdirectionsInnProd
                        * rowColdirectionsInnProd) - (rowRowdirectionsInnProd / 3.0)
                        - LatticeType::discreteVelocityVectors[ALPHA][*rowIndexIncomingVelocity]
                            * LatticeType::discreteVelocityVectors[ALPHA][*rowIndexIncomingVelocity]
                            * (colColdirectionsInnProd - (DIMENSION / 3.0)));

                assert(fabs(kMatrices[contiguousSiteIndex](rowIndex, columnIndex)) < 1e3);
              }

            }
          }

          /**
           * Assemble the L matrix for site contiguousSiteIndex. L is a square matrix (num_incoming_vels x num_incoming_vels)
           *
           * @param contiguousSiteIndex Contiguous site index (for this core)
           */
          inline void AssembleLMatrix(site_t siteLocalIndex)
          {
            unsigned incomingVelsSetSize = incomingVelocities[siteLocalIndex].size();

            lMatrices[siteLocalIndex] = ublas::identity_matrix<distribn_t>(incomingVelsSetSize) + THETA
                * ublas::subrange(kMatrices[siteLocalIndex], 0, incomingVelsSetSize, 0, incomingVelsSetSize);
          }

          /**
           * Compute the LU factorisation of all the L matrices
           */
          inline void FactoriseLMatrices()
          {
            for (std::map<site_t, ublas::matrix<distribn_t> >::iterator iter = lMatrices.begin(); iter
                != lMatrices.end(); ++iter)
            {
              luPermutationMatrices[iter->first]
                  = new ublas::permutation_matrix<std::size_t>(incomingVelocities[iter->first].size());
              int ret = lu_factorize(iter->second, *luPermutationMatrices[iter->first]);
              // If this assertion trips, lMatrices[siteLocalIndex] is singular.
              assert(ret == 0);
            }
          }

          /**
           * Assemble the sigma vector required to assemble the system RHS. We are not including
           * the forcing term used in the paper to drive the flow. This might become necessary
           * for biocolloids.
           *
           * @param contiguousSiteIndex Contiguous site index (for this core)
           * @param sigmaVector sigma vector
           */
          inline void AssembleSigmaVector(const site_t contiguousSiteIndex,
                                          ublas::c_vector<distribn_t, LatticeType::NUMVECTORS>& sigmaVector) // const
          {
            assert(fPostCollision.find(contiguousSiteIndex) != fPostCollision.end());
            assert(fOld.find(contiguousSiteIndex) != fOld.end());
            sigmaVector = fPostCollision[contiguousSiteIndex] - (1 - THETA) * fOld[contiguousSiteIndex];
          }

          /**
           * Assemble the r vector required to assemble the system RHS
           *
           * @param contiguousSiteIndex Contiguous site index (for this core)
           * @param rVector r vector
           */
          inline void AssembleRVector(const site_t contiguousSiteIndex, ublas::vector<distribn_t>& rVector) //const
          {
            ublas::c_vector < distribn_t, LatticeType::NUMVECTORS > sigmaVector;
            AssembleSigmaVector(contiguousSiteIndex, sigmaVector);

            unsigned incomingVelsSetSize = incomingVelocities[contiguousSiteIndex].size();
            unsigned outgoingVelsSetSize = outgoingVelocities[contiguousSiteIndex].size();

            ublas::vector < distribn_t > fNew(outgoingVelsSetSize);

            /*
             *  Assemble a vector with the updated values of the distribution function for the outgoing velocities,
             *  which have already been streamed
             */
            assert(fNew.size() == outgoingVelocities[contiguousSiteIndex].size());
            unsigned index = 0;
            for (std::set<Direction>::const_iterator outgoingDirIter = outgoingVelocities[contiguousSiteIndex].begin(); outgoingDirIter
                != outgoingVelocities[contiguousSiteIndex].end(); ++outgoingDirIter, ++index)
            {
              fNew[index] = *latticeData.GetFNew(contiguousSiteIndex * LatticeType::NUMVECTORS + *outgoingDirIter);
            }

            rVector = THETA * ublas::prod(ublas::subrange(kMatrices[contiguousSiteIndex],
                                                          0,
                                                          incomingVelsSetSize,
                                                          incomingVelsSetSize,
                                                          LatticeType::NUMVECTORS),
                                          fNew) + ublas::prod(kMatrices[contiguousSiteIndex], sigmaVector);
          }
      };
    }
  }
}

#endif /* HEMELB_LB_STREAMERS_JUNKYANG_H */
