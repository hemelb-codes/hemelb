// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#include "redblood/VertexBag.h"
#include "redblood/Interpolation.h"

namespace hemelb
{
  namespace redblood
  {

    namespace detail
    {

      proc_t get_proc(geometry::LatticeData const &latDat, LatticeVector const &pos)
      {
        proc_t procid;
        site_t siteid;
        return latDat.GetContiguousSiteId(pos, procid, siteid) ?
          procid :
          std::numeric_limits<proc_t>::max();
      }

    }

    VertexBag::VertexBag(std::shared_ptr<CellBase const> parent) :
            CellBase(std::shared_ptr<CellBase::CellData>(new CellData(parent->GetTemplateMesh(),
                                                                      parent->GetTag(),
                                                                      parent->GetTemplateName())))
    {
    }
    VertexBag::VertexBag(std::shared_ptr<CellBase const> parent, LatticePosition vertex) :
            CellBase(std::shared_ptr<CellBase::CellData>(new CellData(parent->GetTemplateMesh(),
                                                                      parent->GetTag(),
                                                                      parent->GetTemplateName())))
    {
      addVertex(vertex);
    }

    VertexBag::VertexBag(boost::uuids::uuid const &tag, std::string const &templateName) :
            CellBase(std::make_shared<CellBase::CellData>(Mesh(std::make_shared<MeshData>(MeshData())),
                                                          tag,
                                                          templateName))
    {
    }

    template<class STENCIL>
    std::map<size_t, std::shared_ptr<VertexBag>> splitVertices(
        std::shared_ptr<CellBase const> cell, geometry::LatticeData const &latticeData,
        proc_t selfRegion)
    {
      auto proc_getter = [&latticeData](LatticeVector const &position)
      {
        return detail::get_proc(latticeData, position);
      };
      return splitVertices<STENCIL>(proc_getter, cell, selfRegion);
    }
  }
} // hemelb::redblood
