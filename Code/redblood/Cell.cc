//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//
#include <fstream>
#include <numeric>
#include "redblood/Cell.h"
#include "redblood/Cell.impl.cc"

namespace hemelb
{
  namespace redblood
  {
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
                         moduli.strain,
                         moduli.dilation,
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
                         moduli.strain,
                         moduli.dilation,
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
      for (auto const & neighbors : GetTopology()->facetNeighbors)
      {
        for (auto neighbor : neighbors)
        {
          if (neighbor > static_cast<std::size_t>(current_facet))
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
      std::size_t current_facet(0);
      for (auto const & neighbors : GetTopology()->facetNeighbors)
      {
        for (auto neighbor : neighbors)
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

    std::unique_ptr<CellBase> Cell::cloneImpl() const
    {
      std::unique_ptr<Cell> result(new Cell(GetVertices(),
                                            GetTemplateMesh(),
                                            GetScale(),
                                            GetTemplateName()));
      result->moduli = moduli;
      return std::move(result);
    }

    void writeMesh(std::string const& filename, std::shared_ptr<CellBase const> cell,
                   util::UnitConverter const& converter) {
      log::Logger::Log<log::Debug, log::Singleton>("Writing red blood cell to %s",
          filename.c_str());
      std::ofstream file(filename.c_str());
      assert(file.is_open());
      writeMesh(file, cell, converter);
    }

    void writeMesh(std::ostream & stream, std::shared_ptr<CellBase const> cell,
                   util::UnitConverter const& converter) {
      writeMesh(stream, cell->GetVertices(), cell->GetFacets(), converter);
    }

    void writeVTKMeshWithForces(std::string const &filename, std::shared_ptr<Cell const> cell,
                                util::UnitConverter const &converter)
    {
      log::Logger::Log<log::Debug, log::Singleton>("Writing red blood cell to %s",
                                                   filename.c_str());
      std::ofstream file(filename.c_str());

      if (!file.is_open())
      {
        std::stringstream message;
        message << "Cannot create file '" << filename.c_str() << "', RBC won't be written to disk." << std::endl
                << "Error " << errno << ": " << std::strerror(errno) << std::endl
                << "Bits: " << file.eof() << " " << file.bad() << " " <<  file.fail();
        log::Logger::Log<log::Error, log::OnePerCore>(message.str());
        return;
      }

      writeVTKMeshWithForces(file, cell, converter);
    }

    void writeVTKMeshWithForces(std::ostream &stream, std::shared_ptr<Cell const> cell,
                                util::UnitConverter const &converter)
    {
      PointScalarData point_scalar_data;

      auto push_back_point_data =
          [&point_scalar_data, &cell](const std::string& fieldName, const std::vector<LatticeForceVector>& forces)
          {
            std::vector<LatticeForce> force_magnitudes;
            force_magnitudes.reserve(cell->GetNumberOfNodes());
            for(auto force : forces)
            {
              force_magnitudes.push_back(force.GetMagnitude());
            }

            point_scalar_data.push_back(std::make_pair(fieldName, force_magnitudes));
          };

      auto moduli = cell->moduli;

      {
        std::vector<LatticeForceVector> forces(cell->GetNumberOfNodes(), 0.0);
        cell->facetBending(forces);
        push_back_point_data("bending", forces);
      }

      {
        std::vector<LatticeForceVector> forces(cell->GetNumberOfNodes(), 0.0);
        volumeEnergy(cell->GetVertices(),
                     *cell->GetTemplateMesh().GetData(),
                     moduli.volume,
                     forces,
                     cell->GetScale());
        push_back_point_data("volume", forces);
      }

      {
        std::vector<LatticeForceVector> forces(cell->GetNumberOfNodes(), 0.0);
        surfaceEnergy(cell->GetVertices(),
                      *cell->GetTemplateMesh().GetData(),
                      moduli.surface,
                      forces,
                      cell->GetScale());
        push_back_point_data("surface", forces);
      }

      {
        std::vector<LatticeForceVector> forces(cell->GetNumberOfNodes(), 0.0);
        strainEnergy(cell->GetVertices(),
                     *cell->GetTemplateMesh().GetData(),
                     moduli.strain,
                     moduli.dilation,
                     forces,
                     cell->GetScale());
        push_back_point_data("strain", forces);
      }

      writeVTKMesh(stream,
                   cell->GetVertices(),
                   cell->GetTemplateMesh().GetFacets(),
                   converter,
                   point_scalar_data);
    }

  }
} // hemelb::redblood
