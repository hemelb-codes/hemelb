//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#include <cmath>
#include <iomanip>
#include "redblood/Mesh.h"
#include "redblood/Facet.impl.cc"
#include "constants.h"

// Helper functions in anonymous namespace.
// These are located in separate file so we can easily unit-test them.
namespace hemelb
{
  namespace redblood
  {
    namespace
    {
      // Facet bending energy between two facets
      LatticeEnergy facetBending(Facet const &facetA, Facet const &facetB, Facet const &facetA_eq,
                                  Facet const &facetB_eq, LatticeModulus intensity)
      {
        Angle const theta = orientedAngle(facetA, facetB);
        Angle const theta0 = orientedAngle(facetA_eq, facetB_eq);
        // forces on nodes that are in common
        return 0.5 * intensity * (theta - theta0) * (theta - theta0);
      }

      // Facet bending energy and force between neighboring facets
      LatticeEnergy facetBending(ForceFacet const &facetA, ForceFacet const &facetB,
                                  Facet const &facetA_eq, Facet const &facetB_eq,
                                  LatticeModulus intensity)
      {
        IndexPair const commons = commonNodes(facetA, facetB);
        IndexPair const singles = singleNodes(facetA, facetB);

        LatticePosition const normali = facetA.unitNormal();
        LatticePosition const normalj = facetB.unitNormal();

        // Figures out orientation of the common nodes, whether we have x3 - x1 or x1 - x3.
        const bool orientation =
          facetA(commons.first, commons.second)
            .Cross(facetA(singles.first, commons.second))
            .Dot(normali) < 0e0;

        auto n_ij = normali - normalj * normali.Dot(normalj);
        auto n_ji = normalj - normali * normalj.Dot(normali);
        auto const area_ij = n_ij.GetMagnitude();
        auto const area_ji = n_ji.GetMagnitude();
        if(area_ij > 1e-12)
        {
          n_ij = n_ij / area_ij;
        }
        if(area_ji > 1e-12)
        {
          n_ji = n_ji / area_ji;
        }

        Angle const theta = orientedAngle(facetA, facetB);
        Angle const theta0 = orientedAngle(facetA_eq, facetB_eq);
        const LatticeModulus strength = 5e-1 * intensity * (theta - theta0)
          * (theta < 0e0 ? 1e0: -1e0);
        n_ij = n_ij * (strength/facetB.area());
        n_ji = n_ji * (strength/facetA.area());

        auto const & n1 = facetA(orientation ? commons.first: commons.second);
        auto const & n2 = facetA(singles.first);
        auto const & n3 = facetA(orientation ? commons.second: commons.first);
        auto const & n4 = facetB(singles.second);

        auto & f1 = facetA.GetForce(orientation ? commons.first: commons.second);
        auto & f2 = facetA.GetForce(singles.first);
        auto & f3 = facetA.GetForce(orientation ? commons.second: commons.first);
        auto & f4 = facetB.GetForce(singles.second);

        f1 += (n2 - n3).Cross(n_ji) + (n3 - n4).Cross(n_ij);
        f2 += (n3 - n1).Cross(n_ji);
        f3 += (n1 - n2).Cross(n_ji) + (n4 - n1).Cross(n_ij);
        f4 += (n1 - n3).Cross(n_ij);

        return 0.5 * intensity * (theta - theta0) * (theta - theta0);
      }

      LatticeEnergy facetBending(MeshData::Vertices const &vertices, MeshData const &orig,
                                  size_t facetIndex, size_t neighborIndex, LatticeModulus intensity,
                                  std::vector<LatticeForceVector> &forces)
      {
        return facetBending(ForceFacet(vertices, orig.facets[facetIndex], forces),
                            ForceFacet(vertices, orig.facets[neighborIndex], forces),
                            Facet(orig, facetIndex),
                            Facet(orig, neighborIndex),
                            intensity);
      }
      LatticeEnergy facetBending(MeshData::Vertices const &vertices, MeshData const &orig,
                                  size_t facetIndex, size_t neighborIndex, LatticeModulus intensity)
      {
        return facetBending(Facet(vertices, orig.facets[facetIndex]),
                            Facet(vertices, orig.facets[neighborIndex]),
                            Facet(orig, facetIndex),
                            Facet(orig, neighborIndex),
                            intensity);
      }

      LatticeEnergy volumeEnergy(MeshData::Vertices const &vertices, MeshData const &orig,
                                  LatticeModulus intensity, Dimensionless origMesh_scale = 1e0)
      {
        if (intensity <= 1e-12)
        {
          return 0e0;
        }

        LatticeVolume const vol0 = volume(orig) * origMesh_scale * origMesh_scale * origMesh_scale;
        LatticeVolume const deltaV = volume(vertices, orig.facets) - vol0;
        return intensity * 0.5 * deltaV * deltaV / vol0;
      }

      LatticeEnergy volumeEnergy(MeshData::Vertices const &vertices, MeshData const &orig,
                                  LatticeModulus intensity, std::vector<LatticeForceVector> &forces,
                                  Dimensionless origMesh_scale = 1e0)
      {
        if (intensity <= 1e-12)
        {
          return 0e0;
        }

        assert(orig.vertices.size() == vertices.size());
        LatticeVolume const vol0 = volume(orig) * origMesh_scale * origMesh_scale * origMesh_scale;
        LatticeVolume const deltaV = volume(vertices, orig.facets) - vol0;
        double const strength(intensity / 6.0 * deltaV / vol0);

        for (auto const &facet : orig.facets)
        {
          // Come aliases to make it easier to refer to vertices
          size_t const i0(facet[0]), i1(facet[1]), i2(facet[2]);
          LatticePosition const &a(vertices[i0]);
          LatticePosition const &b(vertices[i1]);
          LatticePosition const &c(vertices[i2]);

          forces[i0] += b.Cross(c) * strength;
          forces[i1] += c.Cross(a) * strength;
          forces[i2] += a.Cross(b) * strength;
        }

        return 0.5 * intensity * deltaV * deltaV / vol0;
      }

      LatticeEnergy surfaceEnergy(MeshData::Vertices const &vertices, MeshData const &orig,
                                   LatticeModulus intensity, Dimensionless origMesh_scale = 1e0)
      {
        LatticeArea const surf0 = area(orig) * origMesh_scale * origMesh_scale;
        LatticeArea const deltaS = area(vertices, orig.facets) - surf0;
        return intensity * 0.5 * deltaS * deltaS / surf0;
      }

      LatticeEnergy surfaceEnergy(MeshData::Vertices const &vertices, MeshData const &orig,
                                   LatticeModulus intensity, std::vector<LatticeForceVector> &forces,
                                   Dimensionless origMesh_scale = 1e0)
      {
        assert(orig.vertices.size() == vertices.size());

        LatticeArea const surf0 = area(orig) * origMesh_scale * origMesh_scale;
        LatticeArea const deltaS = area(vertices, orig.facets) - surf0;
        double const strength = intensity * 0.5 * deltaS / surf0;

        for (size_t facetIndex(0); facetIndex < orig.facets.size(); ++facetIndex)
        {
          ForceFacet facet(vertices, orig.facets, facetIndex, forces);
          LatticePosition const n0 = facet.unitNormal();

          facet.GetForce(0) += n0.Cross(facet(2, 1)) * strength;
          facet.GetForce(1) += n0.Cross(facet(0, 2)) * strength;
          facet.GetForce(2) += n0.Cross(facet(1, 0)) * strength;
        }

        return intensity * 0.5 * deltaS * deltaS / surf0;
      }

      LatticeEnergy strainEnergyDensity(
          std::pair<Dimensionless, Dimensionless> const &strainParams, LatticeModulus shearModulus,
          LatticeModulus dilationModulus)
      {
        Dimensionless const I1 = strainParams.first, I2 = strainParams.second;
        return shearModulus / 12. * (I1 * I1 + 2. * I1 - 2. * I2) + dilationModulus / 12. * I2 * I2;
      }
      LatticeEnergy strainEnergy(Facet const &deformed, Facet const &undeformed,
                                  LatticeModulus shearModulus, LatticeModulus dilationModulus,
                                  Dimensionless origMesh_scale = 1e0)
      {
        return strainEnergyDensity(strainInvariants(deformed, undeformed, origMesh_scale),
                                   shearModulus,
                                   dilationModulus) * undeformed.area();
      }

      LatticeEnergy strainEnergy(ForceFacet const &deformed, Facet const &undeformed,
                                  LatticeModulus shearModulus, LatticeModulus dilationModulus,
                                  Dimensionless origMesh_scale = 1e0)
      {
        // Shape function parameters
        Dimensionless const b0 = undeformed.length(0) * 0.5 * origMesh_scale, b1 =
            (undeformed.length(1) * undeformed.cosine() - undeformed.length(0)) * 0.5
                * origMesh_scale, a1 = -0.5 * undeformed.length(1) * undeformed.sine()
            * origMesh_scale;

        LatticePosition const disps = displacements(deformed, undeformed, origMesh_scale);
        LatticePosition const squaredDisps(squaredDisplacements(disps));
        std::pair<Dimensionless, Dimensionless> const strainInvs(strainInvariants(squaredDisps));
        Dimensionless const I1 = strainInvs.first, I2 = strainInvs.second;
        Dimensionless const w = strainEnergyDensity(strainInvs, shearModulus, dilationModulus);

        // Skalak Parameters
        LatticeModulus const dw_dI1 = shearModulus / 6 * (I1 + 1), dw_dI2 = -shearModulus / 6.
            + dilationModulus / 6. * I2;

        size_t const xx = 0, yy = 1, xy = 2;

        // Derivatives of strain invariants
        Dimensionless const dI1_dGxx = 1., dI1_dGyy = 1., dI2_dGxx = squaredDisps[yy], dI2_dGyy =
            squaredDisps[xx], dI2_dGxy = -2. * squaredDisps[xy];

        // Derivatives of squared deformation tensor
        Dimensionless const dGxx_du1x = 2. * a1 * disps[xx], dGxy_du0x = b0 * disps[xx], dGxy_du1x =
            a1 * disps[xy] + b1 * disps[xx], dGxy_du1y = a1 * disps[yy], dGyy_du0x = 2. * b0
            * disps[xy], dGyy_du0y = 2. * b0 * disps[yy], dGyy_du1x = 2. * b1 * disps[xy],
            dGyy_du1y = 2. * b1 * disps[yy];

        LatticeModulus const force0x = (dw_dI1 * dI1_dGyy * dGyy_du0x
            + dw_dI2 * (dI2_dGyy * dGyy_du0x + dI2_dGxy * dGxy_du0x)), force0y = (dw_dI1 * dI1_dGyy
            * dGyy_du0y + dw_dI2 * dI2_dGyy * dGyy_du0y), force1x = (dw_dI1
            * (dI1_dGxx * dGxx_du1x + dI1_dGyy * dGyy_du1x)
            + dw_dI2 * (dI2_dGxx * dGxx_du1x + dI2_dGyy * dGyy_du1x + dI2_dGxy * dGxy_du1x)),
            force1y = (dw_dI1 * dI1_dGyy * dGyy_du1y
                + dw_dI2 * (dI2_dGyy * dGyy_du1y + dI2_dGxy * dGxy_du1y));

        /// Coordinate system
        LatticePosition const ex = deformed.edge(0).GetNormalised(), ez =
            deformed.edge(0).Cross(deformed.edge(1)).GetNormalised(), ey = ez.Cross(ex);

        LatticeForceVector const force0 = ex * force0x + ey * force0y, force1 = ex * force1x
            + ey * force1y;
        deformed.GetForce(0) -= force0;
        deformed.GetForce(1) -= force1;
        deformed.GetForce(2) += force0 + force1;

        return w * undeformed.area() * origMesh_scale * origMesh_scale;
      }

      LatticeEnergy strainEnergy(MeshData::Vertices const &vertices, MeshData const &origin,
                                 LatticeModulus shearModulus, LatticeModulus dilationModulus,
                                 Dimensionless origMesh_scale = 1e0)
      {
        LatticeEnergy result(0);

        for (size_t i(0); i < origin.facets.size(); ++i)
          result += strainEnergy(Facet(vertices, origin.facets[i]),
                                 Facet(origin, i),
                                 shearModulus,
                                 dilationModulus,
                                 origMesh_scale);

        return result;
      }
      LatticeEnergy strainEnergy(MeshData::Vertices const &vertices, MeshData const &origin,
                                 LatticeModulus shearModulus, LatticeModulus dilationModulus,
                                 std::vector<LatticeForceVector> &forces,
                                 Dimensionless origMesh_scale = 1e0)
      {
        LatticeEnergy result(0);

        for (size_t i(0); i < origin.facets.size(); ++i)
          result += strainEnergy(ForceFacet(vertices, origin.facets[i], forces),
                                 Facet(origin, i),
                                 shearModulus,
                                 dilationModulus,
                                 origMesh_scale);

        return result;
      }
    }
  }
}
