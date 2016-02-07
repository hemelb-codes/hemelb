//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_REDBLOOD_MESH_H
#define HEMELB_REDBLOOD_MESH_H

#include <memory>
#include <array>
#include <vector>
#include <set>
#include <string>
#include <type_traits>
#include <map>
#include "util/Vector3D.h"
#include "util/Matrix3D.h"
#include "util/UnitConverter.h"
#include "units.h"

namespace hemelb
{
  namespace redblood
  {
    //! Holds raw mesh data
    //! Data is separated into vertices and triangular facets
    class MeshData
    {
      public:
        //! Type of containers over indices
        typedef std::array<size_t, 3> Facet;
        //! Facet container type
        typedef std::vector<Facet> Facets;
        //! Vertex container type
        typedef std::vector<LatticePosition> Vertices;
        //! Vertex container
        Vertices vertices;
        //! Facet container
        Facets facets;
    };
    static_assert(std::is_standard_layout<MeshData>::value, "Needed for MPI data type");
    static_assert(
        std::is_default_constructible<MeshData>::value
        and std::is_nothrow_default_constructible<MeshData>::value
        and std::is_move_constructible<MeshData>::value
        and std::is_nothrow_move_constructible<MeshData>::value
        and std::is_copy_constructible<MeshData>::value
        and std::is_copy_assignable<MeshData>::value
        and (not std::is_nothrow_copy_assignable<MeshData>::value)
        and (not std::is_pod<MeshData>::value),
        "Explicit type characteristics"
    );

    LatticePosition barycenter(MeshData const &mesh);
    LatticePosition barycenter(MeshData::Vertices const &vertices);
    LatticeVolume volume(MeshData const &mesh);
    LatticeVolume volume(MeshData::Vertices const &vertices, MeshData::Facets const &facets);
    LatticeArea area(MeshData const &mesh);
    LatticeArea area(MeshData::Vertices const &vertices, MeshData::Facets const &facets);
    //! Orients facet inward, or inward
    void orientFacets(MeshData &mesh, bool outward = true);

    //! Holds raw topology data
    class MeshTopology
    {
      public:
        //! Type for map from vertices to facets
        typedef std::vector<std::set<std::size_t> > VertexToFacets;
        //! Type for map from facets to its neighbors
        typedef std::vector<std::array<std::size_t, 3> > FacetNeighbors;
        //! For each vertex, lists the facet indices
        VertexToFacets vertexToFacets;
        //! For each facet, lists the neighboring facets
        FacetNeighbors facetNeighbors;

        // Creates mesh topology from mesh data
        MeshTopology(MeshData const &mesh);
        // Empty topology
        MeshTopology()
        {
        }
    };
    static_assert(
        std::is_default_constructible<MeshTopology>::value
        and (not std::is_nothrow_default_constructible<MeshTopology>::value)
        and std::is_move_constructible<MeshTopology>::value
        and std::is_nothrow_move_constructible<MeshTopology>::value
        and std::is_copy_constructible<MeshTopology>::value
        and std::is_copy_assignable<MeshTopology>::value
        and (not std::is_nothrow_copy_assignable<MeshTopology>::value)
        and std::is_standard_layout<MeshTopology>::value
        and (not std::is_pod<MeshTopology>::value),
        "Explicit type characteristics"
    );

    //! Triangular mesh
    class Mesh
    {
        //! Differentiates between two shallow and deep copy constructors
        struct deepcopy_tag
        {
        };

        //! Deep copy constructor
        Mesh(Mesh const &c, deepcopy_tag const &) :
            mesh(new MeshData(*c.mesh)), topology(new MeshTopology(*c.topology))
        {
        }

      public:
        //! Initializes mesh from mesh data
        Mesh(std::shared_ptr<MeshData> const &mesh) :
            mesh(mesh), topology(new MeshTopology(*mesh))
        {
        }
        //! Initializes mesh from mesh data and topology
        Mesh(std::shared_ptr<MeshData> const &mesh, std::shared_ptr<MeshTopology> const &topo) :
            mesh(mesh), topology(topo)
        {
        }
        //! Initialize mesh by copying data
        Mesh(MeshData const &data) :
            mesh(new MeshData(data)), topology(new MeshTopology(data))
        {
        }
        //! Shallow copy constructor
        Mesh(Mesh const &in) :
            mesh(in.mesh), topology(in.topology)
        {
        }

        //! Determines barycenter of mesh
        LatticePosition GetBarycenter() const
        {
          return barycenter(*mesh);
        }
        //! Computes volume of the mesh
        LatticeVolume GetVolume() const
        {
          return volume(*mesh);
        }
        //! Computes area of the mesh
        LatticeArea GetArea() const
        {
          return area(*mesh);
        }
        //! Connectivity data
        std::shared_ptr<const MeshTopology> GetTopology() const
        {
          return topology;
        }
        //! Returns mesh data
        std::shared_ptr<MeshData> GetData()
        {
          return mesh;
        }
        //! Returns mesh data
        std::shared_ptr<MeshData const> GetData() const
        {
          return mesh;
        }
        //! Returns mesh data
        void SetData(MeshData const &data)
        {
          mesh.reset(new MeshData(data));
        }

        //! Makes a copy of the mesh
        Mesh clone() const
        {
          return Mesh(*this, deepcopy_tag());
        }

        //! Scale mesh around barycenter
        void operator*=(Dimensionless const &scale);
        //! Scale by matrix around barycenter
        void operator*=(util::Matrix3D const &rotation);
        //! Translate mesh
        void operator+=(LatticePosition const &offset);
        //! Transform mesh
        void operator+=(std::vector<LatticePosition> const &displacements);

        //! Number of nodes
        size_t GetNumberOfNodes() const
        {
          return mesh->vertices.size();
        }
        //! Const access to vertices
        MeshData::Vertices &GetVertices()
        {
          return mesh->vertices;
        }
        //! Const access to vertices
        MeshData::Vertices const &GetVertices() const
        {
          return mesh->vertices;
        }
        //! Const access to vertices
        MeshData::Vertices::const_reference GetVertex(size_t site) const
        {
          return mesh->vertices[site];
        }
        //! Const access to facets
        MeshData::Facets const &GetFacets() const
        {
          return mesh->facets;
        }

        //! True if both mesh share the same data
        bool isSameData(Mesh const &other) const
        {
          return mesh == other.mesh and topology == other.topology;
        }

      protected:
        //! Holds actual data about the mesh
        std::shared_ptr<MeshData> mesh;
        //! Holds topology information
        std::shared_ptr<MeshTopology> topology;
    };
    static_assert(
        (not std::is_default_constructible<Mesh>::value)
        and (not std::is_nothrow_default_constructible<Mesh>::value)
        and std::is_move_constructible<Mesh>::value
        and (not std::is_nothrow_move_constructible<Mesh>::value)
        and std::is_copy_constructible<Mesh>::value
        and std::is_copy_assignable<Mesh>::value
        and std::is_nothrow_copy_assignable<Mesh>::value
        and std::is_standard_layout<Mesh>::value
        and (not std::is_pod<Mesh>::value),
        "Explicit type characteristics"
    );

    //! Read mesh from file
    //! Format is from T. Krueger's thesis
    std::shared_ptr<MeshData> readMesh(std::string const &filename);
    //! Read mesh from file
    //! Format is from T. Krueger's thesis
    std::shared_ptr<MeshData> readMesh(std::istream &stream);
    //! Write mesh from file
    //! Format is from T. Krueger's thesis
    void writeMesh(std::ostream &stream, MeshData const &data, util::UnitConverter const &);
    //! Write mesh from file
    //! Format is from T. Krueger's thesis
    void writeMesh(std::string const &filename, MeshData const &data, util::UnitConverter const &);
    //! Write mesh to file in VTK XML format
    void writeVTKMesh(std::ostream &, MeshData const &, util::UnitConverter const &);
    //! Write mesh to file in VTK XML format
    void writeVTKMesh(std::string const &, MeshData const &, util::UnitConverter const&);
    //! Write mesh to file in VTK XML format
    void writeVTKMesh(std::ostream &, MeshData::Vertices const &, MeshData::Facets const &,
                      util::UnitConverter const&);

    //! Tetrahedron of a depth
    //! Depth refers to the number of triangular subdivision in each facet
    Mesh tetrahedron(unsigned int depth = 0);
    //! Flat triangular shape
    //! Depth refers to the number of triangular subdivision in each facet
    //! There are two initial facets, referring to the same three vertices.
    Mesh pancakeSamosa(unsigned int depth = 0);
    //! Refine a mesh by decomposing each facet into four triangles
    Mesh refine(Mesh mesh, unsigned int depth = 0);
    //! Creates an ico sphere (regular triangular mesh over a sphere)
    Mesh icoSphere(unsigned int depth = 0);

    //! Orients facet outward, or inward
    void orientFacets(Mesh &mesh, bool outward = true);
  }
}
#endif
