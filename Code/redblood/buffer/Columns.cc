//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#include <iterator>
#include "redblood/buffer/Columns.h"
#include "util/Matrix3D.h"

namespace hemelb
{
  namespace redblood
  {
    namespace buffer
    {
      namespace
      {
        //! Rotation matrix for a to b
        //! Rotation matrix maps one basis onto another:
        //!    - a0 (normalised a) maps to b0 (normalised b)
        //!    - a0.Cross(b0) maps to a0.Cross(b0)
        //!    - a0.Cross(a0, Cross(b0)) maps to b0.Cross(a0.Cross(b0))
        util::Matrix3D rotMat(LatticePosition const& a, LatticePosition const& b)
        {
          util::Matrix3D result;
          LatticePosition const a0 = a.GetNormalised();
          LatticePosition const b0 = b.GetNormalised();
          // Special case where a0 == b0: Rotation is identity
          // Also works for a = 0 and b = 0... But that's just bad input
          if(a0.Cross(b0).GetMagnitude() < 1e-8)
          {
            for(size_t i(0); i < 3; ++i)
            {
              for(size_t j(0); j < 3; ++j)
              {
                result[i][j] = i == j ? 1: 0;
              }
            }
            return result;
          }

          LatticePosition const u = (a0.Cross(b0)).GetNormalised();
          LatticePosition const a1 = a0.Cross(u).GetNormalised();
          LatticePosition const b1 = b0.Cross(u).GetNormalised();

          // HemeLB doesn't need linear algebra, so it doesn't have it - cos that would be
          // complicated - which in turn implies it doesn't need it. Instead, do lets use a square
          // wheel.
          util::Matrix3D A, B;
          A[0][0] = a0[0]; A[0][1] = a0[1]; A[0][2] = a0[2];
          A[1][0] =  u[0]; A[1][1] =  u[1]; A[1][2] =  u[2];
          A[2][0] = a1[0]; A[2][1] = a1[1]; A[2][2] = a1[2];
          B[0][0] = b0[0]; B[0][1] =  u[0]; B[0][2] = b1[0];
          B[1][0] = b0[1]; B[1][1] =  u[1]; B[1][2] = b1[1];
          B[2][0] = b0[2]; B[2][1] =  u[2]; B[2][2] = b1[2];
          return B * A;
        }

        LatticeDistance maxExtension(
            MeshData::Vertices const &vertices, LatticePosition const &direction)
        {
          auto const barycenter = hemelb::redblood::barycenter(vertices);
          LatticeDistance result(0);
          auto const normalised = direction.GetNormalised();
          auto maxdist = [&normalised, &result, &barycenter](LatticePosition const &b)
          {
            result = std::max(result, std::abs(normalised.Dot(b - barycenter)));
          };
          std::for_each(vertices.begin(), vertices.end(), maxdist);
          return 2e0 * result;
        }

        LatticePosition maxExtensions(
            MeshData::Vertices const &vertices,
            LatticePosition const &col,
            LatticePosition const& normal)
        {
          LatticeDistance const x = maxExtension(vertices, normal);
          LatticeDistance const y = maxExtension(vertices, col.Cross(normal));
          LatticeDistance const z = maxExtension(vertices, col);
          return LatticePosition(x, y, z);
        }
      }

#     ifndef HEMELB_DOING_UNITTESTS
      ColumnPositionIterator::ColumnPositionIterator(
              std::shared_ptr<Cylinder const> cylinder, MeshData::Vertices const& vertices,
              LatticePosition cellAxis, LatticePosition colAxis,
              LatticeDistance separation)
          : IndexIterator(LatticeVector(0, 0, 0), LatticeVector(0, 0, 0)), cylinder(cylinder)
      {
        if(cellAxis.GetMagnitude() < 1e-8)
        {
          throw Exception() << "Cell axis cannot be 0";
        }
        if(colAxis.GetMagnitude() < 1e-8)
        {
          throw Exception() << "Column axis cannot be 0";
        }
        if(std::abs(colAxis.Dot(cylinder->normal)) > 1e-8) 
        {
          throw Exception() << "Column axis should be perpendicular to cylinder axis";
        }

        // Rotation is opposite to the one that will be applied to the mesh
        auto const antiRot = rotMat(colAxis, cellAxis);
        auto const extents
          = maxExtensions(vertices, antiRot * colAxis, antiRot * cylinder->normal) + separation;

        max.x = std::numeric_limits<LatticeCoordinate>::max();
        max.y = static_cast<LatticeCoordinate>(std::ceil(cylinder->radius / extents.y));
        max.z = static_cast<LatticeCoordinate>(std::ceil(cylinder->radius / extents.z));
        min.x = 0; min.y = -max.y; min.z = -max.z;
        current = min;

        major = colAxis.GetNormalised() * extents.z;
        minor = colAxis.GetNormalised().Cross(cylinder->normal.GetNormalised()) * extents.y;
        depth = cylinder->normal.GetNormalised() * extents.x;
        spacing = (extents + separation) / 2.0;

        if(not IsInside()) operator++();
      }

      bool ColumnPositionIterator::IsInside() const
      {
        // avoids missing defs for operator^
        PhysicalDistance const a
          =  major.GetMagnitudeSquared() / std::pow(cylinder->radius - spacing.z, 2e0);
        PhysicalDistance const b
          = minor.GetMagnitudeSquared() / std::pow(cylinder->radius - spacing.y, 2e0);
        return (current.z * current.z) * a + (current.y * current.y) * b < 1e0;
      }

      void ColumnPositionIterator::operator++()
      {
        assert(IsValid());
        IndexIterator::operator++();

#       ifndef NDEBUG
        size_t breaker = 100, iter = 0;
#       endif
        while(not IsInside())
        {
          IndexIterator::operator++();
#         ifndef NDEBUG
            if((++iter) >= breaker)
            {
              throw Exception() << "Probably in an infinite loop";
            }
#         endif
        }
      }

      ColumnCellDrop::ColumnCellDrop(
          std::shared_ptr<Cylinder const> cylinder, CellContainer::value_type cell,
          LatticePosition cellAxis, LatticePosition colAxis,
          LatticeDistance wallSeparation)
        : iterator(cylinder, cell->GetVertices(), cellAxis, colAxis, wallSeparation),
          templateCell(cell->clone())
      {
        // Cell is rotated to correct orientation
        *templateCell *= rotMat(cellAxis, colAxis);
        // And centered at zero
        *templateCell -= templateCell->GetBarycenter();
      }

      CellContainer::value_type ColumnCellDrop::operator()()
      {
        ++iterator;
        std::shared_ptr<CellBase> result(std::move(templateCell->clone()));
        *result += *iterator;
        return result;
      }


#     endif
    } // namespace buffer
  }
} // hemelb::redblood
