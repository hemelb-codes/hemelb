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
#include "redblood/Particle.h"
#include "unittests/redblood/Fixtures.h"

namespace hemelb { namespace unittests {

class ParticleTests : public EnergyVsGradientFixture {
    CPPUNIT_TEST_SUITE(ParticleTests);
    CPPUNIT_TEST(testParticleEnergy);
    CPPUNIT_TEST_SUITE_END();

    struct ParticleEnergy {
      mutable redblood::Particle particle;
      ParticleEnergy(redblood::Mesh const &_mesh,
          redblood::Mesh const _template)
        : particle(_mesh, _template) {
          particle.moduli.bending = 0.888;
          particle.moduli.surface = 1.127;
          particle.moduli.volume = 1.015231;
          particle.moduli.strain = 1.047524;
          particle.moduli.dilation = 0.945524;
      }
      PhysicalEnergy operator()(redblood::MeshData const& _mesh) const {
        particle.SetData(redblood::MeshData(_mesh));
        return particle();
      }
      PhysicalEnergy operator()(
          redblood::MeshData const &_mesh,
          std::vector<LatticeForceVector> &_forces) const {
        particle.SetData(redblood::MeshData(_mesh));
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

    void testParticleEnergy() {
      energyVsForces(ParticleEnergy(mesh, original));
    }

protected:
    redblood::MeshData original;
};

CPPUNIT_TEST_SUITE_REGISTRATION(ParticleTests);
}}

#endif // ONCE

