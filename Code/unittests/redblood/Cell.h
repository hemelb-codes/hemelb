//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_REDBLOOD_PARTICLE_H
#define HEMELB_UNITTESTS_REDBLOOD_PARTICLE_H

#include <cppunit/TestFixture.h>
#include "redblood/Cell.h"
#include "unittests/redblood/Fixtures.h"

namespace hemelb { namespace unittests { namespace redblood {

class CellTests : public EnergyVsGradientFixture {
    CPPUNIT_TEST_SUITE(CellTests);
    CPPUNIT_TEST(testCellEnergy);
    CPPUNIT_TEST_SUITE_END();

    struct CellEnergy {
      mutable Cell particle;
      CellEnergy(Mesh const &_mesh, Mesh const _template)
        : particle(_mesh, _template) {
          particle.moduli.bending = 0.888;
          particle.moduli.surface = 1.127;
          particle.moduli.volume = 1.015231;
          particle.moduli.strain = 1.047524;
          particle.moduli.dilation = 0.945524;
      }
      PhysicalEnergy operator()(MeshData const& _mesh) const {
        particle.SetData(MeshData(_mesh));
        return particle();
      }
      PhysicalEnergy operator()(
          MeshData const &_mesh,
          std::vector<LatticeForceVector> &_forces) const {
        particle.SetData(MeshData(_mesh));
        return particle(_forces);
      }
    };
public:

    void setUp() {
      BasisFixture::setUp();
      original = mesh;
      mesh.vertices[0] += LatticePosition(-0.01, 0.02342, 0.03564);
      mesh.vertices[1] += LatticePosition(0.0837, -0.012632, 0.0872935);
      mesh.vertices[2] += LatticePosition(0.02631, -0.00824223, -0.098362);
    }

    void testCellEnergy() {
      energyVsForces(CellEnergy(mesh, original));
    }

protected:
    MeshData original;
};

CPPUNIT_TEST_SUITE_REGISTRATION(CellTests);
}}}

#endif // ONCE

