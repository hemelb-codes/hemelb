// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_REDBLOOD_CELLTESTS_H
#define HEMELB_UNITTESTS_REDBLOOD_CELLTESTS_H

#include <cppunit/TestFixture.h>
#include "redblood/Cell.h"
#include "unittests/redblood/Fixtures.h"

namespace hemelb
{
  namespace unittests
  {
    namespace redblood
    {
      class CellTests : public EnergyVsGradientFixture
      {
          CPPUNIT_TEST_SUITE (CellTests);
          CPPUNIT_TEST (testCellCopyShallowness);
          CPPUNIT_TEST (testCellClone);
          CPPUNIT_TEST (testCellEnergy);
          CPPUNIT_TEST (testNullTemplateScaling);
          CPPUNIT_TEST (testTemplateScaling);
          CPPUNIT_TEST (testGetAverageEdgeLength);CPPUNIT_TEST_SUITE_END();

          struct CellEnergy
          {
              mutable Cell particle;
              CellEnergy(Mesh const &mesh, Mesh const origMesh) :
                  particle(mesh, origMesh, 2.5) // Set scale to != 1 to test area scaling in energy calculations
              {
                particle.moduli.bending = 0.888;
                particle.moduli.surface = 1.127;
                particle.moduli.volume = 1.015231;
                particle.moduli.strain = 1.047524;
                particle.moduli.dilation = 0.945524;
              }
              LatticeEnergy operator()(MeshData const &mesh) const
              {
                particle.GetVertices() = mesh.vertices;
                return particle();
              }
              LatticeEnergy operator()(MeshData const &mesh,
                                       std::vector<LatticeForceVector> &forces) const
              {
                particle.GetVertices() = mesh.vertices;
                return particle(forces);
              }
          };

	template <typename VEC3>
	auto vector_of_zero_vec3(std::size_t count) {
	  return std::vector<VEC3>{count, VEC3::Zero()};
	}

        public:
          void setUp()
          {
            BasisFixture::setUp();
            original = mesh;
            mesh.vertices[0] += LatticePosition(-0.01, 0.02342, 0.03564);
            mesh.vertices[1] += LatticePosition(0.0837, -0.012632, 0.0872935);
            mesh.vertices[2] += LatticePosition(0.02631, -0.00824223, -0.098362);
          }

          void testCellEnergy()
          {
            energyVsForces(CellEnergy(mesh, original));
          }

          void testNullTemplateScaling()
          {
            Mesh templateMesh(original);
	    auto forces = vector_of_zero_vec3<LatticeForceVector>(original.vertices.size());

            std::vector<Dimensionless> scales;
            scales.push_back(1.0);
            scales.push_back(0.8);
            scales.push_back(1.2);
            for (auto const scale : scales)
            {
              auto other = templateMesh.clone();
              other *= scale;
              auto cell = GetCellWithEnergy(other, templateMesh, scale);

              CPPUNIT_ASSERT_DOUBLES_EQUAL(cell.GetScale(), scale, 1e-8);
              auto const zero = cell(forces);
              CPPUNIT_ASSERT_DOUBLES_EQUAL(zero, 0e0, 1e-8);

              for (auto const &force : forces)
              {
                CPPUNIT_ASSERT(helpers::is_zero(force));
              }
            }
          }

          void testTemplateScaling()
          {
            Mesh templateMesh(original);
            auto scaled = templateMesh.clone();
            scaled.GetVertices()[0] += LatticePosition(-0.01, 0.02342, 0.03564);
            scaled.GetVertices()[1] += LatticePosition(0.0837, -0.012632, 0.0872935);
            scaled.GetVertices()[2] += LatticePosition(0.02631, -0.00824223, -0.098362);

            auto uforces = vector_of_zero_vec3<LatticeForceVector>(original.vertices.size());
	    auto sforces = vector_of_zero_vec3<LatticeForceVector>(original.vertices.size());

            scaled *= 1.1;
            auto const uenergy = GetCellWithEnergy(scaled, templateMesh, 1.1)(uforces);
            auto scaledTemplateMesh = templateMesh.clone();
            scaledTemplateMesh *= 1.1;
            auto const senergy = GetCellWithEnergy(scaled, scaledTemplateMesh)(sforces);

            CPPUNIT_ASSERT_DOUBLES_EQUAL(uenergy, senergy, 1e-12);
            auto i_unscaled = uforces.cbegin();

            for (auto scaled_force : sforces)
            {
              CPPUNIT_ASSERT(helpers::is_zero(* (i_unscaled++) - scaled_force));
            }
          }

          void testCellCopyShallowness()
          {
            Cell cell0(original);
            Cell cell1(cell0);
            CPPUNIT_ASSERT(cell0.GetTemplateMesh().GetData() == cell1.GetTemplateMesh().GetData());
            CPPUNIT_ASSERT(&cell0.GetVertices() != &cell1.GetVertices());
            CPPUNIT_ASSERT(cell0.GetTag() != cell1.GetTag());

            Cell cell2(cell0, CellBase::shallow_clone());
            CPPUNIT_ASSERT(cell0.GetTemplateMesh().GetData() == cell2.GetTemplateMesh().GetData());
            CPPUNIT_ASSERT(&cell0.GetVertices() == &cell2.GetVertices());
            CPPUNIT_ASSERT(cell0.GetTag() == cell2.GetTag());
          }

          void testCellClone()
          {
            Cell cell0(original);
            auto cell1 = cell0.clone();
            CPPUNIT_ASSERT(cell0.GetTemplateMesh().GetData() == cell1->GetTemplateMesh().GetData());
            CPPUNIT_ASSERT(&cell0.GetVertices() != &cell1->GetVertices());
            CPPUNIT_ASSERT(cell0.GetTag() != cell1->GetTag());
          }

          void testGetAverageEdgeLength()
          {
            // MeshData original is a tetrahedron with 3 edges of length 1 and 3 edges of length sqrt(2)
            Cell cell(original);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(cell.GetAverageEdgeLength(), (1 + sqrt(2)) / 2, 1e-12);
          }

        protected:
          MeshData original;

          Cell GetCellWithEnergy(Mesh const &a, Mesh const &b, Dimensionless s = 1e0) const
          {
            Cell cell(a, b, s);
            cell.moduli.bending = 0.888;
            cell.moduli.surface = 1.127;
            cell.moduli.volume = 1.015231;
            cell.moduli.strain = 1.047524;
            cell.moduli.dilation = 0.945524;
            return cell;
          }
      };
      CPPUNIT_TEST_SUITE_REGISTRATION (CellTests);
    }
  }
}

#endif  // ONCE
