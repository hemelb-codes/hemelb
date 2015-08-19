//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//
#include <fstream>
#include "redblood/Cell.h"
// Helper functions in anonymous namespace.
// These are located in separate file so we can easily unit-test them.
#include "redblood/Cell.impl.cc"

// Contains definition of CellBase::CellData
// Should not be included by anyone except here and VertexBag
// In fact, there are no header guards on this file, just to make sure it is not included wrongly.
#include "redblood/CellDataDefinition.impl.cc"

namespace hemelb
{
  namespace redblood
  {
    CellBase::CellBase(MeshData::Vertices &&verticesIn, Mesh const &origMesh, Dimensionless scaleIn,
                       std::string const & templateName) :
        data(new CellData(std::move(verticesIn), origMesh, scaleIn, templateName))
    {
    }
    CellBase::CellBase(MeshData::Vertices const &verticesIn, Mesh const &origMesh,
                       Dimensionless scaleIn, std::string const & templateName) :
        data(new CellData(verticesIn, origMesh, scaleIn, templateName))
    {
    }
    CellBase::CellBase(CellBase const& cell) :
        data(new CellData(*cell.data))
    {
    }
    CellBase::CellBase(Mesh const &mesh, Mesh const &origMesh, Dimensionless scaleIn,
                       std::string const & templateName) :
        data(new CellData(mesh.GetVertices(), origMesh, scaleIn, templateName))
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
    std::string const & CellBase::GetTemplateName() const
    {
      return data->templateName;
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

    LatticeEnergy Cell::operator()() const
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
    LatticeEnergy Cell::operator()(std::vector<LatticeForceVector> &forces) const
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

    LatticeEnergy Cell::facetBending() const
    {
      if (std::abs(moduli.bending) < 1e-8)
      {
        return 0e0;
      }

      LatticeEnergy result(0);
      site_t current_facet(0);
      for (auto const & neighbors: GetTopology()->facetNeighbors)
      {
        for (auto neighbor: neighbors)
        {
          if (neighbor > current_facet)
          {
            result += hemelb::redblood::facetBending(data->vertices,
                                                     *data->templateMesh.GetData(),
                                                     current_facet,
                                                     neighbor,
                                                     moduli.bending);
          }
        }
        ++current_facet;
      }

      return result;
    }

    LatticeEnergy Cell::facetBending(std::vector<LatticeForceVector> &forces) const
    {
      if (std::abs(moduli.bending) < 1e-8)
      {
        return 0e0;
      }

      LatticeEnergy result(0);
      typedef MeshTopology::FacetNeighbors::const_iterator FacetIterator;
      site_t current_facet(0);
      for (auto const & neighbors: GetTopology()->facetNeighbors)
      {
        for (auto neighbor: neighbors)
        {
          if (neighbor > current_facet)
          {
            result += hemelb::redblood::facetBending(data->vertices,
                                                     *data->templateMesh.GetData(),
                                                     current_facet,
                                                     neighbor,
                                                     moduli.bending,
                                                     forces);
          }
        }
        ++current_facet;
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
      std::unique_ptr<Cell> result(new Cell(GetVertices(),
                                            GetTemplateMesh(),
                                            GetScale(),
                                            GetTemplateName()));
      result->moduli = moduli;
      return std::move(result);
    }

    boost::uuids::uuid const & CellBase::GetTag() const
    {
      return data->tag;
    }

    void writeVTKMesh(
        std::string const &filename, std::shared_ptr<CellBase const> cell,
        util::UnitConverter const &converter)
    {
      log::Logger::Log<log::Debug, log::Singleton>("Writing red blood cell to %s",
                                                   filename.c_str());
      std::ofstream file(filename.c_str());
      writeVTKMesh(file, cell, converter);
    }

    void writeVTKMesh(
        std::ostream &stream, std::shared_ptr<CellBase const> cell,
        util::UnitConverter const &converter)
    {
      writeVTKMesh(stream, cell->GetVertices(), cell->GetTemplateMesh().GetFacets(), converter);
    }

#   ifndef NDEBUG
    void checkCellDataCharacteristics()
    {
      static_assert(
          (not std::is_default_constructible<CellBase::CellData>::value)
          and (not std::is_nothrow_default_constructible<CellBase::CellData>::value)
          and std::is_move_constructible<CellBase::CellData>::value
          and (not std::is_nothrow_move_constructible<CellBase::CellData>::value)
          and std::is_copy_constructible<CellBase::CellData>::value
          and (not std::is_copy_assignable<CellBase::CellData>::value)
          and (not std::is_nothrow_copy_assignable<CellBase::CellData>::value)
          and std::is_standard_layout<CellBase::CellData>::value
          and (not std::is_pod<CellBase::CellData>::value),
          "Explicit type characteristics"
      );
    }
#   endif
  }
} // hemelb::redblood
