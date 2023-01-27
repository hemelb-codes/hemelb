// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <catch2/catch.hpp>

#include "redblood/Cell.h"
#include "redblood/Mesh.h"
#include "redblood/MeshIO.h"
#include "resources/Resource.h"
#include "tests/redblood/Fixtures.h"

namespace hemelb
{
  namespace tests
  {
    redblood::KruegerMeshIO io = {};
    class SadCellTests
    {
    protected:
      redblood::Cell cell;
      
    public:
      SadCellTests() :
	cell(
	     io.readFile(resources::Resource("sad.msh").Path(), true)->vertices,
	     io.readFile(resources::Resource("rbc_ico_1280.msh").Path(), true)
	     )
      {
	// sad mesh is in physical units, move to something with fewer decimals.
	cell *= 1e0 / 4.1e-6;
      }

#     define HEMELB_MACRO(name)			\
      void test_ ## name()			\
      {						\
	cell.moduli.name = 1e0;		\
	checkNumericalForces();			\
	cell.moduli.name = 0e0;		\
      }						\
      void test_ ## name ## ForceDirection()	\
      {						\
	cell.moduli.name = 1e0;		\
	checkForceDirection();			\
	cell.moduli.name = 0e0;		\
      }

      HEMELB_MACRO(bending);
      HEMELB_MACRO(surface);
      HEMELB_MACRO(volume);
      HEMELB_MACRO(dilation);
      HEMELB_MACRO(strain);
#     undef HEMELB_MACRO

      void checkForceDirection()
      {
	std::vector<LatticeForceVector> forces(cell.GetVertices().size(), LatticeForceVector::Zero());
	auto const e0 = cell(forces);
	for (size_t i(0); i < forces.size(); ++i) {
	  cell.GetVertices()[i] += forces[i] * 1e-6;
	}
	auto const e1 = cell();
	REQUIRE(e0 > e1);
      }

      void checkNumericalForces()
      {
	std::vector<LatticeForceVector> forces(cell.GetVertices().size(), LatticeForceVector::Zero());
	cell(forces);
	for (size_t i(0); i < cell.GetVertices().size(); ++i) {
	  for (int j(0); j < 4; ++j) {
	    auto const epsilon = 1e-3;
	    LatticePosition const direction = LatticePosition(
	      j == 0 ? 1 : (j >= 3 ? random() : 0),
	      j == 1 ? 1 : (j >= 3 ? random() : 0),
	      j == 2 ? 1 : (j >= 3 ? random() : 0)
	    ).GetNormalised();

	    double const expected = derivative(direction, i, epsilon);
	    double const actual = -Dot(forces[i], direction);
	    REQUIRE(Approx(expected).margin(1e-7).epsilon(1e-4) == actual);
	  }
	}
      }

      double derivative(LatticePosition const &dir, size_t i, double epsilon)
      {
	double result = 0e0;
	auto const orig = cell.GetVertices()[i];
	auto &vertex = cell.GetVertices()[i];
	std::vector<double> const coeffs { -1. / 280., 4. / 105., -1. / 5., 4. / 5. };
	for (size_t j(0); j < coeffs.size(); ++j) {
	  vertex = orig + dir * double(coeffs.size() - j) * epsilon;
	  result += coeffs[j] * cell();
	  vertex = orig - dir * double(coeffs.size() - j) * epsilon;
	  result -= coeffs[j] * cell();
	}
	vertex = orig;
	return result / epsilon;
      }

    };

#   define HEMELB_MACRO(name)						\
    METHOD_AS_TEST_CASE(SadCellTests::test_ ## name,			\
			"Sad cell test " #name,				\
			"[redblood][.long]");				\
    METHOD_AS_TEST_CASE(SadCellTests::test_ ## name ## ForceDirection,	\
			"Sad cell test " # name " force direction",	\
     			"[redblood][.long]");

    HEMELB_MACRO(bending);
    HEMELB_MACRO(surface);
    HEMELB_MACRO(volume);
    HEMELB_MACRO(dilation);
    HEMELB_MACRO(strain);
#   undef HEMELB_MACRO
  }
}

