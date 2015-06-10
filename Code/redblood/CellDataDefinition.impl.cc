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

namespace hemelb
{
  namespace redblood
  {
    // Advanced declaration
    class VertexBag;
    // Holds the data of CellBase
    // Although its own data members are public, it can only be accessed by CellBase and derived
    // classes. So encapsulation is not broken. Indeed, the definition of this class should not be
    // available to anyone but the CellBase class family. Adding Get/Setters would mean more code
    // for the exact same functionality.
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

      private:
        // Gives access to special constructor for VertexBag
        friend class VertexBag;
        // Constructor that can copy the tag
        // Should only be used by people in the know.
        // It breaks the unicity of the cell.
        CellData(Mesh const &origMesh, boost::uuids::uuid const &uuid)
          : templateMesh(origMesh), scale(1e0), tag(uuid)
        {
          assert(scale > 1e-12);
        }
    };

#   ifndef NDEBUG
      void check_cell_data_characteristics() {
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
}
