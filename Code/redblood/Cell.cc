//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>

#include "redblood/Cell.h"
// Helper functions in anonymous namespace.
// These are located in separate file so we can easily unit-test them.
#include "redblood/Cell.impl.cc"

namespace hemelb
{
  namespace redblood
  {
    class CellBase::CellData
    {
      public:
        CellData(MeshData::Vertices &&verticesIn, Mesh const &origMesh, Dimensionless scaleIn = 1e0) :
            vertices(std::move(verticesIn)), templateMesh(origMesh), scale(scaleIn),
            tag(boost::uuids::random_generator()())
        {
          assert(scale > 1e-12);
        }
        CellData(MeshData::Vertices const &verticesIn, Mesh const &origMesh, Dimensionless scaleIn =
                     1e0) :
            vertices(verticesIn), templateMesh(origMesh), scale(scaleIn),
            tag(boost::uuids::random_generator()())
        {
          assert(scale > 1e-12);
        }
        CellData(CellData const& c)
          : vertices(c.vertices), templateMesh(c.templateMesh), scale(c.scale),
            tag(boost::uuids::random_generator()())
        {
        }
        CellData(CellData && c)
          : vertices(std::move(c.vertices)), templateMesh(std::move(c.templateMesh)),
            scale(c.scale), tag(std::move(c.tag))
        {
        }
        //! Holds list of vertices for this cell
        MeshData::Vertices vertices;
        //! Unmodified original mesh
        Mesh templateMesh;
        //! Scale factor for the template;
        Dimensionless scale;
        //! Uuid tag
        boost::uuids::uuid const tag;
    };

    CellBase::CellBase(MeshData::Vertices &&verticesIn, Mesh const &origMesh, Dimensionless scaleIn) :
        data(new CellData(std::move(verticesIn), origMesh, scaleIn))
    {
    }
    CellBase::CellBase(MeshData::Vertices const &verticesIn, Mesh const &origMesh,
                       Dimensionless scaleIn) :
        data(new CellData(verticesIn, origMesh, scaleIn))
    {
    }
    CellBase::CellBase(CellBase const& cell) :

        data(new CellData(*cell.data))
    {
    }
    CellBase::CellBase(Mesh const &mesh, Mesh const &origMesh, Dimensionless scaleIn) :
        data(new CellData(mesh.GetVertices(), origMesh, scaleIn))
    {
    }
    CellBase::CellBase(Mesh const &mesh) :
        data(new CellData(mesh.GetVertices(), mesh.clone()))
    {
    }
    CellBase::CellBase(std::shared_ptr<MeshData> const &mesh) :
        data(new CellData(mesh->vertices, Mesh(*mesh)))
    {
    }

    void CellBase::operator=(Mesh const &mesh)
    {
      data->templateMesh = mesh;
      data->vertices = data->templateMesh.GetVertices();
      data->scale = 1e0;
    }

    //! Unmodified mesh
    Mesh const &CellBase::GetTemplateMesh() const
    {
      return data->templateMesh;
    }
    //! Facets for the mesh
    MeshData::Facets const &CellBase::GetFacets() const
    {
      return data->templateMesh.GetData()->facets;
    }
    //! Vertices of the cell
    MeshData::Vertices const &CellBase::GetVertices() const
    {
      return data->vertices;
    }
    //! Vertices of the cell
    MeshData::Vertices &CellBase::GetVertices()
    {
      return data->vertices;
    }
    //! Topology of the (template) mesh
    std::shared_ptr<MeshTopology const> CellBase::GetTopology() const
    {
      return data->templateMesh.GetTopology();
    }
    site_t CellBase::GetNumberOfNodes() const
    {
      return static_cast<site_t>(data->vertices.size());
    }
    MeshData::Vertices::value_type CellBase::GetBarycenter() const
    {
      return barycenter(data->vertices);
    }

    //! Scale to apply to the template mesh
    void CellBase::SetScale(Dimensionless scaleIn)
    {
      assert(data->scale > 1e-12);
      data->scale = scaleIn;
    }
    //! Scale to apply to the template mesh
    Dimensionless CellBase::GetScale() const
    {
      return data->scale;
    }

    PhysicalEnergy Cell::operator()() const
    {
      return facetBending() // facet bending unaffected by template scale
      + volumeEnergy(data->vertices, *data->templateMesh.GetData(), moduli.volume, data->scale)
          + surfaceEnergy(data->vertices,
                          *data->templateMesh.GetData(),
                          moduli.surface,
                          data->scale)
          + strainEnergy(data->vertices,
                         *data->templateMesh.GetData(),
                         moduli.dilation,
                         moduli.strain,
                         data->scale);
    }
    PhysicalEnergy Cell::operator()(std::vector<LatticeForceVector> &forces) const
    {
      assert(forces.size() == data->vertices.size());
      return facetBending(forces)
          + volumeEnergy(data->vertices,
                         *data->templateMesh.GetData(),
                         moduli.volume,
                         forces,
                         data->scale)
          + surfaceEnergy(data->vertices,
                          *data->templateMesh.GetData(),
                          moduli.surface,
                          forces,
                          data->scale)
          + strainEnergy(data->vertices,
                         *data->templateMesh.GetData(),
                         moduli.dilation,
                         moduli.strain,
                         forces,
                         data->scale);
    }

    PhysicalEnergy Cell::facetBending() const
    {
      if (std::abs(moduli.bending) < 1e-8)
      {
        return 0e0;
      }

      PhysicalEnergy result(0);
      typedef MeshTopology::FacetNeighbors::const_iterator FacetIterator;
      FacetIterator i_facet = GetTopology()->facetNeighbors.begin();
      FacetIterator const i_facetEnd = GetTopology()->facetNeighbors.end();

      for (size_t current(0); i_facet != i_facetEnd; ++i_facet, ++current)
      {
        for (size_t i(0); i < 3; ++i)
          if ( (*i_facet)[i] > current)
            result += hemelb::redblood::facetBending(data->vertices,
                                                     *data->templateMesh.GetData(),
                                                     current,
                                                     (*i_facet)[i],
                                                     moduli.bending);
      }

      return result;
    }

    PhysicalEnergy Cell::facetBending(std::vector<LatticeForceVector> &forces) const
    {
      if (std::abs(moduli.bending) < 1e-8)
      {
        return 0e0;
      }

      PhysicalEnergy result(0);
      typedef MeshTopology::FacetNeighbors::const_iterator FacetIterator;
      FacetIterator i_facet = GetTopology()->facetNeighbors.begin();
      FacetIterator const i_facetEnd = GetTopology()->facetNeighbors.end();

      for (size_t current(0); i_facet != i_facetEnd; ++i_facet, ++current)
      {
        for (size_t i(0); i < 3; ++i)
          if ( (*i_facet)[i] > current)
            result += hemelb::redblood::facetBending(data->vertices,
                                                     *data->templateMesh.GetData(),
                                                     current,
                                                     (*i_facet)[i],
                                                     moduli.bending,
                                                     forces);
      }

      return result;
    }

    void CellBase::operator*=(Dimensionless const &scaleIn)
    {
      auto const barycenter = GetBarycenter();

      for (auto &vertex : data->vertices)
      {
        vertex = (vertex - barycenter) * scaleIn + barycenter;
      }
    }
    void CellBase::operator*=(util::Matrix3D const &rotation)
    {
      auto const barycenter = GetBarycenter();

      for (auto &vertex : data->vertices)
      {
        rotation.timesVector(vertex - barycenter, vertex);
        vertex += barycenter;
      }
    }
    void CellBase::operator+=(LatticePosition const &offset)
    {
      for (auto &vertex : data->vertices)
      {
        vertex += offset;
      }
    }
    void CellBase::operator+=(std::vector<LatticePosition> const &displacements)
    {
      assert(displacements.size() == data->vertices.size());
      auto i_disp = displacements.begin();

      for (auto &vertex : data->vertices)
      {
        vertex += * (i_disp++);
      }
    }

    std::unique_ptr<CellBase> Cell::cloneImpl() const
    {
      std::unique_ptr<Cell> result(new Cell(GetVertices(), GetTemplateMesh(), GetScale()));
      result->moduli = moduli;
      return std::move(result);
    }

    boost::uuids::uuid const & CellBase::GetTag() const
    {
      return data->tag;
    }

  }
} // hemelb::redblood
