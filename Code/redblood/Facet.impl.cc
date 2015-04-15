#ifndef HEMELB_REDBLOOD_FACET_H
#define HEMELB_REDBLOOD_FACET_H

#include <cmath>
#include <iomanip>
#include "constants.h"
#include "redblood/Mesh.h"

// Helper functions in anonymous namespace.
// These are located in separate file so we can easily unit-test them.
// They are in an anonymous namespace since they likely not needed as a public
// API.
namespace hemelb
{
  namespace redblood
  {
    namespace
    {
      // Helper class to avoid explicit indexing over vertices
      class Facet
      {
        public:
          // References nodes of a facet
          LatticePosition const *nodes[3];
          // Indices of nodes in original array
          MeshData::Facet const &indices;
          Facet(MeshData const &mesh, size_t index) :
              Facet(mesh.vertices, mesh.facets[index])
          {
          }
          Facet(MeshData::Vertices const &vertices, MeshData::Facets const &facets, size_t index) :
              Facet(vertices, facets[index])
          {
          }
          Facet(MeshData::Vertices const &vertices, MeshData::Facet const &indices) :
              indices(indices)
          {
            nodes[0] = &vertices[indices[0]];
            nodes[1] = &vertices[indices[1]];
            nodes[2] = &vertices[indices[2]];
          }

          // returns an edge nodes[i] - nodes[j]
          LatticePosition operator()(size_t i, size_t j) const;
          // returns node i
          LatticePosition const &operator()(size_t i) const;
          // Returns edges
          LatticePosition edge(size_t i) const;
          // Returns edge length
          PhysicalDistance length(size_t i) const;
          // Returns angle cosine
          Dimensionless cosine() const;
          // Returns angle sine
          Dimensionless sine() const;
          // Computes vector normal to facet
          // Order of edges determines direction of normal
          // This is taken straight from T. Krueger's code
          LatticePosition normal() const;
          // Unit vector normal to facet
          LatticePosition unitNormal() const;
          // Area of the facet
          PhysicalArea area() const;
      };

      // Facet that also includes forces
      class ForceFacet : public Facet
      {
        public:
          // References forces on a node
          LatticeForceVector *forces[3];
          ForceFacet(MeshData::Vertices const &vertices, MeshData::Facet const &indices,
                     std::vector<LatticeForceVector> &forcesIn) :
              Facet(vertices, indices)
          {
            forces[0] = &forcesIn[indices[0]];
            forces[1] = &forcesIn[indices[1]];
            forces[2] = &forcesIn[indices[2]];
          }
          ForceFacet(MeshData const &mesh, size_t index, std::vector<LatticeForceVector> &forcesIn) :
              ForceFacet(mesh.vertices, mesh.facets[index], forcesIn)
          {
          }
          ForceFacet(MeshData::Vertices const &vertices, MeshData::Facets const &facets,
                     size_t index, std::vector<LatticeForceVector> &forcesIn) :
              ForceFacet(vertices, facets[index], forcesIn)
          {
          }
          LatticeForceVector &GetForce(size_t i) const
          {
            return * (forces[i]);
          }
      };

      LatticePosition Facet::operator()(size_t i, size_t j) const
      {
        return (*this)(i) - (*this)(j);
      }
      LatticePosition const &Facet::operator()(size_t i) const
      {
        return * (nodes[i]);
      }
      LatticePosition Facet::edge(size_t i) const
      {
        switch (i)
        {
          case 0:
            return operator()(2, 1);

          case 1:
            return operator()(0, 1);

          case 2:
            return operator()(2, 0);

          default:
            return LatticePosition(0, 0, 0);
        };
      }
      PhysicalDistance Facet::length(size_t i) const
      {
        return edge(i).GetMagnitude();
      }
      Dimensionless Facet::cosine() const
      {
        return edge(0).Dot(edge(1)) / (length(0) * length(1));
      }
      Dimensionless Facet::sine() const
      {
        // 0 < angle < 180 in a triangle, so sine always positive in case of
        // interest. I think.
        return edge(0).Cross(edge(1)).GetMagnitude() / (length(0) * length(1));
      }
      LatticePosition Facet::normal() const
      {
        LatticePosition const edgeA(operator()(0, 1));
        LatticePosition const edgeB(operator()(2, 1));
        return edgeA.Cross(edgeB);
      }
      LatticePosition Facet::unitNormal() const
      {
        return normal().Normalise();
      }
      PhysicalArea Facet::area() const
      {
        return normal().GetMagnitude() * 0.5;
      }

      // Just a helper function to check whether v is in a
      bool contains(MeshData::Facet const &a, MeshData::Facet::value_type v)
      {
        return a[0] == v or a[1] == v or a[2] == v;
      }

      // Computes common nodes for neighboring facets
      // This routine will report nonsense if facets are not neighbors
      typedef std::pair<size_t, size_t> IndexPair;
      IndexPair commonNodes(Facet const &a, Facet const &b)
      {
        // First node differs, hence other two nodes in common
        if (not contains(b.indices, a.indices[0]))
        {
          return IndexPair(1, 2);
        }

        // First node in common, second node differs
        if (not contains(b.indices, a.indices[1]))
        {
          return IndexPair(0, 2);
        }

        // First node and second in common
        return IndexPair(0, 1);
      }
      // Returns common edge with specific direction
      LatticePosition commonEdge(Facet const &a, Facet const &b)
      {
        IndexPair common(commonNodes(a, b));
        return a(common.first, common.second);
      }

      // Figures out nodes that are not in common
      // Returns non-sense if the nodes are not neighbors.
      IndexPair singleNodes(Facet const &a, Facet const &b)
      {
        IndexPair result;
        size_t mappingB[3] = { 4, 4, 4 };

        for (size_t i(0); i < 3; ++i)
          if (a.indices[i] == b.indices[0])
          {
            mappingB[0] = i;
          }
          else if (a.indices[i] == b.indices[1])
          {
            mappingB[1] = i;
          }
          else if (a.indices[i] == b.indices[2])
          {
            mappingB[2] = i;
          }
          else
          {
            result.first = i;
          }

        for (size_t i(0); i < 3; ++i)
          if (mappingB[i] == 4)
          {
            result.second = i;
            break;
          }

        return result;
      }

      // Computes angle between two facets
      Angle angle(LatticePosition const &a, LatticePosition const &b)
      {
        Angle const cosine(a.Dot(b));

        if (cosine >= (1e0 - 1e-6))
        {
          return 0e0;
        }
        else if (cosine <= - (1e0 - 1e-6))
        {
          return PI;
        }

        return std::acos(cosine);
      }
      Angle angle(Facet const &a, Facet const &b)
      {
        return angle(a.unitNormal(), b.unitNormal());
      }

      // Angle with orientation
      // Computes angle between two vectors, including orientation.
      // The orientation should be a vector with an out-of-plane component (eg
      // parallel to the cross product of the two normals)
      Angle orientedAngle(LatticePosition const &a, LatticePosition const &b,
                          LatticePosition const &orient)
      {
        Angle const result(angle(a, b));
        return a.Cross(b).Dot(orient) <= 0e0 ?
          result :
          -result;
      }
      Angle orientedAngle(Facet const &a, Facet const &b, LatticePosition const &orient)
      {
        return orientedAngle(a.unitNormal(), b.unitNormal(), orient);
      }
      Angle orientedAngle(Facet const &a, Facet const &b)
      {
        return orientedAngle(a, b, commonEdge(a, b));
      }

      // Returns Dxx, Dyy, Dxy packed in vector
      LatticePosition displacements(Facet const &deformed, Facet const &ref,
                                    Dimensionless origMesh_scale = 1e0)
      {
        PhysicalDistance const dlength0 = deformed.length(0), rlength0 = ref.length(0)
            * origMesh_scale, dlength1 = deformed.length(1), rlength1 = ref.length(1)
            * origMesh_scale;
        Dimensionless const dsine = deformed.sine(), rsine = ref.sine();
        return LatticePosition(
        // Dxx
                               dlength0 / rlength0,
                               // Dyy
                               (dlength1 * dsine) / (rlength1 * rsine),
                               // Dxy
                               (dlength1 / rlength1 * deformed.cosine()
                                   - dlength0 / rlength0 * ref.cosine()) / rsine);
      }

      // Returns Gxx, Gyy, Gxy packed in vector
      LatticePosition squaredDisplacements(LatticePosition const &disp)
      {
        return LatticePosition(
        // Gxx
                               disp[0] * disp[0],
                               // Gyy
                               disp[2] * disp[2] + disp[1] * disp[1],
                               // Gxy
                               disp[0] * disp[2]);
      }
      LatticePosition squaredDisplacements(Facet const &deformed, Facet const &ref,
                                           Dimensionless origMesh_scale = 1e0)
      {
        return squaredDisplacements(displacements(deformed, ref, origMesh_scale));
      }

      // Strain invariants I1 and I2
      std::pair<Dimensionless, Dimensionless> strainInvariants(LatticePosition const &squaredDisp)
      {
        return std::pair<Dimensionless, Dimensionless>(squaredDisp[0] + squaredDisp[1] - 2.0,
                                                       squaredDisp[0] * squaredDisp[1]
                                                           - squaredDisp[2] * squaredDisp[2] - 1e0);
      }

      std::pair<Dimensionless, Dimensionless> strainInvariants(Facet const &deformed,
                                                               Facet const &ref,
                                                               Dimensionless origMesh_scale = 1e0)
      {
        return strainInvariants(squaredDisplacements(deformed, ref, origMesh_scale));
      }
    }
  }
}

#endif
