// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <fstream>
#include <numeric>
#include "redblood/Cell.h"
#include "redblood/CellEnergy.h"

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

  }
} // hemelb::redblood
