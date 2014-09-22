//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_REDBLOOD_FIXTURES_H
#define HEMELB_UNITTESTS_REDBLOOD_FIXTURES_H
namespace hemelb { namespace unittests {


// Helper functions for tests.
struct Comparisons {
  static bool is_zero(util::Vector3D<double> const &_in, double _tol = 1e-8) {
    return std::sqrt(_in.GetMagnitudeSquared()) < _tol;
  }
  static bool is_zero(double const _in, double _tol = 1e-8) {
    return std::abs(_in) < _tol;
  }
};

class TetrahedronFixture : public CppUnit::TestFixture, public Comparisons {
  public:
    void setUp() {
      // facets at something degrees from one another
      mesh.vertices.push_back(LatticePosition(0, 0, 0));
      mesh.vertices.push_back(LatticePosition(1, 0, 0));
      mesh.vertices.push_back(LatticePosition(0, 1, 0));
      mesh.vertices.push_back(LatticePosition(0, 0, 1));


      redblood::MeshData::t_Facet indices;
      indices[0] = 0; indices[1] = 1; indices[2] = 2;
      mesh.facets.push_back(indices);
      indices[0] = 0; indices[1] = 2; indices[2] = 3;
      mesh.facets.push_back(indices);
      indices[0] = 0; indices[1] = 3; indices[2] = 1;
      mesh.facets.push_back(indices);
      indices[0] = 1; indices[1] = 3; indices[2] = 2;
      mesh.facets.push_back(indices);
    }

  protected:
    hemelb::redblood::MeshData mesh;
};

}}
#endif
