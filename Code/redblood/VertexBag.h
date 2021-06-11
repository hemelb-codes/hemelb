// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_REDBLOOD_VERTEXBAG_H
#define HEMELB_REDBLOOD_VERTEXBAG_H

#include <map>
#include <string>
#include <limits>
#include <boost/uuid/uuid.hpp>

#include "redblood/Cell.h"
#include "redblood/Interpolation.h"
#include "redblood/stencil.h"
#include "geometry/LatticeData.h"
#include "Exception.h"

namespace hemelb
{
  namespace redblood
  {

    //! \brief Holds vertices from another cell
    //! \details However it mimics CellBase so that it can be used in circumstances similar to other
    //! CellBase objects.
    class VertexBag : public CellBase
    {
      public:
        //! Creates a bag to refer to vertices from parent
        VertexBag(std::shared_ptr<CellBase const> parent);
        //! Convenience constructor that also adds first vertex
        VertexBag(std::shared_ptr<CellBase const> parent, LatticePosition vertex);
        //! Construct from scratch
        VertexBag(boost::uuids::uuid const &tag, std::string const &templateName);
        //! Adds a vertex
        //! Shorter than going through GetVertices.
        void addVertex(LatticePosition const &pos)
        {
          CellBase::GetVertices().push_back(pos);
        }

        virtual LatticeEnergy operator()() const override
        {
          throw Exception() << "This object is a fake.";
        }

        virtual LatticeEnergy operator()(std::vector<LatticeForceVector> &in) const override
        {
          throw Exception() << "This object is a fake.";
        }
        std::unique_ptr<VertexBag> clone() const
        {
          return std::unique_ptr<VertexBag>(static_cast<VertexBag*>(cloneImpl().release()));
        }
      private:
        std::unique_ptr<CellBase> cloneImpl() const override
        {
          throw Exception() << "This object cannot be cloned";
          return std::unique_ptr<CellBase>();
        }
    };

    //! Splits vertices from cell depending on which region it belongs to
    //! \return map where the key refers to the region, and the value holds the vertices for that
    //!         region.
    //! \param[in] cell for which to split vertices
    //! \param[in] latticeData holds information about geometry and regions over which to split
    //! \param[in] selfRegion does not store information about this region. Generally, will refer to
    //!             current processor. By default, a value so large, that if it refers to an actual
    //!             region, you've got other problems.
    std::map<size_t, std::shared_ptr<VertexBag>> splitVertices(
        std::shared_ptr<CellBase const> cell, geometry::LatticeData const &latticeData,
        proc_t selfRegion = std::numeric_limits<proc_t>::max(), redblood::stencil::types stencil =
            redblood::stencil::types::FOUR_POINT);

    namespace detail {
      proc_t get_proc(geometry::LatticeData const &latDat, LatticeVector const &pos);

      //! Set of procs affected by this position
      //! \param[in] latDat will tell us which site belongs to which proc
      //! \param[in] iterator a  stencil iterator going over affected lattice points
      template<class STENCIL, class T_FUNC>
      std::set<proc_t> procsAffectedByPosition(
          T_FUNC get_proc, InterpolationIterator<STENCIL> &&iterator, proc_t avoid
					       ) {
        std::set<proc_t> result;
        for (; iterator.IsValid(); ++iterator)
        {
          auto const procid = get_proc(*iterator);
          if (procid != avoid)
          {
            result.insert(procid);
          }
        }
        return result;
      }

      //! Set of procs affected by this position
      //! \param[in] latDat will tell us which site belongs to which proc
      //! \param[in] position for which to figure out affected processes
      //! \param[in] stencil giving interaction range
      template<class STENCIL, class T_FUNC>
      std::set<proc_t> procsAffectedByPosition(
          T_FUNC get_proc, LatticePosition const &position,
          proc_t avoid = std::numeric_limits<proc_t>::max()
					       ) {
        return procsAffectedByPosition(get_proc, interpolationIterator<STENCIL>(position), avoid);
      }

      template<class STENCIL, class T_FUNC>
      std::map<size_t, std::shared_ptr<VertexBag>> splitVertices(
          T_FUNC get_proc, std::shared_ptr<CellBase const> cell,
          proc_t avoid = std::numeric_limits<proc_t>::max())
      {
        std::map<size_t, std::shared_ptr<VertexBag>> result;
        for (auto const &vertex : cell->GetVertices())
        {
          auto const regions = procsAffectedByPosition<STENCIL>(get_proc, vertex, avoid);
          for (auto const region : regions)
          {
            auto const i_bag = result.find(region);
            if (i_bag == result.end())
            {
              result.emplace(region, std::make_shared<VertexBag>(cell, vertex));
            }
            else
            {
              i_bag->second->addVertex(vertex);
            }
          }
        }
        return result;
      }
    }
  }
} // namespace hemelb::redblood
#endif
