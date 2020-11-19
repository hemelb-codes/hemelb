// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <cstring>
#include <sstream>

#include "lb/kernels/Kernels.h"
#include "lb/kernels/rheologyModels/RheologyModels.h"
#include "lb/kernels/momentBasis/DHumieresD3Q15MRTBasis.h"
#include "lb/kernels/momentBasis/DHumieresD3Q19MRTBasis.h"

#include "tests/lb/LbTestsHelper.h"

#include "tests/helpers/FourCubeLatticeData.h"
#include "tests/helpers/FourCubeBasedTestFixture.h"

namespace hemelb
{
  namespace tests
  {
    // File contains tests for the functionality of the
    // lattice-Boltzmann kernels. This includes the functions for
    // calculating hydrodynamic variables and for performing
    // collisions.
    template <typename T>
    struct CollisionTester : public helpers::FourCubeBasedTestFixture {
      using LATTICE = lb::lattices::D3Q15;
      static constexpr auto NV = LATTICE::NUMVECTORS;
      using KERNEL = T;
      using HYDRO = lb::kernels::HydroVars<KERNEL>;
      using DISTS = distribn_t[NV];
      using VEC = util::Vector3D<distribn_t>;
      const distribn_t allowedError = 1e-10;
      KERNEL kernel;
      DISTS f_original;
      CollisionTester() : kernel(initParams) {
	// Initialise the original f distribution to something asymmetric.
	LbTestsHelper::InitialiseAnisotropicTestData<LATTICE>(0, f_original);
      }
      void make_feq(const distribn_t& density, const VEC& momentum, DISTS& feq) const;
      void make_fpostcol(const distribn_t* feq, DISTS& fpc) const;
      void do_tests(HYDRO& hydroVars, const distribn_t& expectedDensity, const VEC& expectedMomentum) {
	DISTS expectedFEq;
	make_feq(expectedDensity, expectedMomentum, expectedFEq);

	// Now compare the expected and actual values
	LbTestsHelper::CompareHydros(expectedDensity,
				     expectedMomentum,
				     expectedFEq,
				     hydroVars,
				     allowedError);
	// Do the collision and test the result.
	kernel.DoCollide(lbmParams, hydroVars);

	// Get the expected post-collision densities.
	DISTS expectedPostCollision;
	make_fpostcol(hydroVars.GetFEq().f, expectedPostCollision);

	// Compare.
	for (unsigned int ii = 0; ii < NV; ++ii) {
	  REQUIRE(Approx(expectedPostCollision[ii]).margin(allowedError) == hydroVars.GetFPostCollision()[ii]);
	}
      }
    };
    template<>
    void CollisionTester<lb::kernels::EntropicAnsumali<lb::lattices::D3Q15>>::make_feq(const distribn_t& density, const VEC& momentum, DISTS& feq) const {
      LbTestsHelper::CalculateAnsumaliEntropicEqmF<LATTICE>(density,
							    momentum,
							    feq);
    }
    template<>
    void CollisionTester<lb::kernels::EntropicAnsumali<lb::lattices::D3Q15>>::make_fpostcol(const distribn_t* feq, DISTS& fpc) const {
      LbTestsHelper::CalculateEntropicCollision<LATTICE>(f_original,
							 feq,
							 lbmParams->GetTau(),
							 lbmParams->GetBeta(),
							 fpc);
    }

    template<>
    void CollisionTester<lb::kernels::EntropicChik<lb::lattices::D3Q15>>::make_feq(const distribn_t& density, const VEC& momentum, DISTS& feq) const {
      LATTICE::CalculateEntropicFeqChik(density,
					momentum.x, momentum.y, momentum.z,
					feq);
    }
    template<>
    void CollisionTester<lb::kernels::EntropicChik<lb::lattices::D3Q15>>::make_fpostcol(const distribn_t* feq, DISTS& fpc) const {
      LbTestsHelper::CalculateEntropicCollision<LATTICE>(f_original,
							 feq,
							 lbmParams->GetTau(),
							 lbmParams->GetBeta(),
							 fpc);
    }

    template<>
    void CollisionTester<lb::kernels::LBGK<lb::lattices::D3Q15>>::make_feq(const distribn_t& density, const VEC& momentum, DISTS& feq) const {
      LbTestsHelper::CalculateLBGKEqmF<lb::lattices::D3Q15>(density,
							    momentum,
							    feq);
    }
    template<>
    void CollisionTester<lb::kernels::LBGK<lb::lattices::D3Q15>>::make_fpostcol(const distribn_t* feq, DISTS& fpc) const {
      LbTestsHelper::CalculateLBGKCollision<lb::lattices::D3Q15>(f_original,
								 feq,
								 lbmParams->GetOmega(),
								 fpc);
    }

    TEMPLATE_TEST_CASE_METHOD(CollisionTester, "KernelTests - Ansumali & Chikatamarla entropic, LBGK calculations and collision", "[lb][kernels]",
			      lb::kernels::EntropicAnsumali<lb::lattices::D3Q15>,
			      lb::kernels::EntropicChik<lb::lattices::D3Q15>,
			      lb::kernels::LBGK<lb::lattices::D3Q15>) {
      using Fix = CollisionTester<TestType>;

      SECTION("use the function that calculates density, momentum and f_eq") {
	typename Fix::HYDRO hydroVars(this->f_original);
      	// Calculate density, momentum, equilibrium f.
      	this->kernel.CalculateDensityMomentumFeq(hydroVars, 0);

      	// Calculate expected values
      	distribn_t expectedDensity = 12.0; // (sum 1 to 15) / 10
      	auto expectedMomentum = LbTestsHelper::CalculateMomentum<typename Fix::LATTICE>(hydroVars.f);
	this->do_tests(hydroVars, expectedDensity, expectedMomentum);
      }

      SECTION("use the function that leaves density and momentum and  calculates f_eq") {
      	typename Fix::HYDRO hydroVars(this->f_original);
      	// Manually set density and momentum and calculate eqm f.
      	hydroVars.density = 1.0;
      	hydroVars.momentum = typename Fix::VEC{0.4, 0.5, 0.6};
      	this->kernel.CalculateFeq(hydroVars, 1);

      	// Calculate expected values
      	distribn_t expectedDensity = 1.0; // Should be unchanged
      	typename Fix::VEC expectedMomentum{ 0.4, 0.5, 0.6 };
      	this->do_tests(hydroVars, expectedDensity, expectedMomentum);
      }
    }

    TEST_CASE_METHOD(helpers::FourCubeBasedTestFixture, "LBGKNNCalculationsAndCollision") {
      using LATTICE = lb::lattices::D3Q15;
      static constexpr auto NV = LATTICE::NUMVECTORS;
      using RHEO_MODEL = lb::kernels::rheologyModels::CarreauYasudaRheologyModelHumanFit;
      using KERNEL = lb::kernels::LBGKNN<RHEO_MODEL, LATTICE>;
      using HYDRO = lb::kernels::HydroVars<KERNEL>;
      using DISTS = distribn_t[NV];
      using VEC = util::Vector3D<distribn_t>;
      const distribn_t numTolerance = 1e-10;
      auto apprx = [&](double x) {
	return Approx(x).margin(numTolerance);
      };

      // We need two kernel instances if we want to work with two
      // different sets of data (and keep the computed values of tau
      // consistent). One to be used with CalculateDensityMomentumFeq
      // and another with CalculateFeq.
      KERNEL lbgknn0(initParams), lbgknn1(initParams);

      // When testing this streamer, it is important to consider that
      // tau is defined per site.  Use two different sets of initial
      // conditions across the domain to check that different
      // shear-rates and relaxation times are computed and stored
      // properly.
      // 
      // Using {f_, velocities}setA for odd site indices and {f_,
      // velocities}setB for the even ones
      // A => 0, B=>1

      DISTS f_original[2];
      // distribn_t* f_original;
      for (unsigned int ii = 0; ii < NV; ++ii) {
	f_original[0][ii] = (1 + ii) / 10.0;
	f_original[1][ii] = (NV - ii) / 10.0;
      }

      HYDRO hydroVars0[2] = {f_original[0], f_original[1]};
      HYDRO hydroVars1[2] = {f_original[0], f_original[1]};

      VEC momentum[2] = {VEC{ 0.4, 0.5, 0.6 }, VEC{ -0.4, -0.5, -0.6 }};

      for (site_t site_index = 0; site_index < numSites; site_index++) {
	// Test part 1: Equilibrium function, density, and momentum
	// are computed identically to the standard LBGK. Local
	// relaxation times are implicitly computed by
	// CalculateDensityMomentumFeq

	// Case 0: test the kernel function for calculating density,
	// momentum and f_eq.
	// 
	// Case 1: test the function that uses a given density and
	// momentum, and calculates f_eq.
	const unsigned set = site_index % 2;

	// Calculate density, momentum, equilibrium f.
	lbgknn0.CalculateDensityMomentumFeq(hydroVars0[set], site_index);

	// Manually set density and momentum and calculate eqm f.
	hydroVars1[set].density = 1.0;
	hydroVars1[set].momentum = momentum[set];

	lbgknn1.CalculateFeq(hydroVars1[set], site_index);

	// Calculate expected values.
	distribn_t expectedDensity0 = 12.0; // (sum 1 to 15) / 10
	distribn_t expectedDensity1 = 1.0; // Unchanged

	auto expectedMomentum0 = LbTestsHelper::CalculateMomentum<LATTICE>(hydroVars0[set].f);
	auto&& expectedMomentum1 = momentum[set];

	DISTS expectedFEq0;
	LbTestsHelper::CalculateLBGKEqmF<LATTICE>(expectedDensity0,
						  expectedMomentum0,
						  expectedFEq0);
	DISTS expectedFEq1;
	LbTestsHelper::CalculateLBGKEqmF<lb::lattices::D3Q15>(expectedDensity1,
							      expectedMomentum1,
							      expectedFEq1);

	// Now compare the expected and actual values.
	LbTestsHelper::CompareHydros(expectedDensity0,
				     expectedMomentum0,
				     expectedFEq0,
				     hydroVars0[set],
				     numTolerance);
	LbTestsHelper::CompareHydros(expectedDensity1,
				     expectedMomentum1,
				     expectedFEq1,
				     hydroVars1[set],
				     numTolerance);

	// Test part 2: Test that the array containing the local
	// relaxation times has the right length and test against some
	// hardcoded values.  Correctness of the relaxation time
	// calculator is tested in RheologyModelTest.h

	// A second call to the Calculate* functions will make sure
	// that the newly computed tau is used in DoCollide as
	// opposite to the default Newtonian tau used during the first
	// time step.
	lbgknn0.CalculateDensityMomentumFeq(hydroVars0[set], site_index);
	lbgknn1.CalculateFeq(hydroVars1[set], site_index);

	distribn_t computedTau0 = hydroVars0[set].tau;
	REQUIRE(numSites == (site_t) lbgknn0.GetTauValues().size());

	distribn_t expectedTau0 = set
	  ? 0.50009102385
	  : 0.50009217276;
	REQUIRE(apprx(expectedTau0) == computedTau0);

	distribn_t computedTau1 = hydroVars1[set].tau;
	REQUIRE(numSites == (site_t) lbgknn1.GetTauValues().size());

	distribn_t expectedTau1 = set
	  ? 0.50009010316
	  : 0.50009016145;

	REQUIRE(apprx(expectedTau1) == computedTau1);

	// Test part 3: Collision depends on the local relaxation time
	// 
	// Do the collision and test the result.
	lbgknn0.DoCollide(lbmParams, hydroVars0[set]);
	lbgknn1.DoCollide(lbmParams, hydroVars1[set]);

	// Get the expected post-collision densities.
	DISTS expectedPostCollision0;
	DISTS expectedPostCollision1;

	distribn_t localOmega0 = -1.0 / computedTau0;
	distribn_t localOmega1 = -1.0 / computedTau1;

	LbTestsHelper::CalculateLBGKCollision<LATTICE>(f_original[set],
						       hydroVars0[set].GetFEq().f,
						       localOmega0,
						       expectedPostCollision0);

	LbTestsHelper::CalculateLBGKCollision<LATTICE>(f_original[set],
						       hydroVars1[set].GetFEq().f,
						       localOmega1,
						       expectedPostCollision1);

	// Compare.
	for (unsigned int ii = 0; ii < lb::lattices::D3Q15::NUMVECTORS; ++ii) {
	  REQUIRE(apprx(expectedPostCollision0[ii]) == hydroVars0[set].GetFPostCollision()[ii]);
	  REQUIRE(apprx(expectedPostCollision1[ii]) == hydroVars1[set].GetFPostCollision()[ii]);
	}
      }
    }

    template <typename L, typename B>
    struct MRTTestFixture : public helpers::FourCubeBasedTestFixture {
      using LATTICE = L;
      static constexpr auto NV = LATTICE::NUMVECTORS;
      using BASIS = B;
      using KERNEL = lb::kernels::MRT<BASIS>;
      using HYDRO = lb::kernels::HydroVars<KERNEL>;
      using DISTS = distribn_t[NV];
      using VEC = util::Vector3D<distribn_t>;

      MRTTestFixture() {
	const distribn_t allowedError = 1e-10;
	KERNEL mrtLbgkEquivalentKernel(initParams);

	// Simulate LBGK by relaxing all the MRT modes to equilibrium
	// with the same time constant.
	std::vector<distribn_t> relaxationParameters;
	distribn_t oneOverTau = 1.0 / lbmParams->GetTau();
	relaxationParameters.resize(BASIS::NUM_KINETIC_MOMENTS, oneOverTau);
	mrtLbgkEquivalentKernel.SetMrtRelaxationParameters(relaxationParameters);

	// Initialise the original f distribution to something asymmetric.
	DISTS f_original;
	LbTestsHelper::InitialiseAnisotropicTestData<LATTICE>(0, f_original);
	HYDRO hydroVars0(f_original);

	// Calculate density, momentum, equilibrium f.
	mrtLbgkEquivalentKernel.CalculateDensityMomentumFeq(hydroVars0, 0);

	// Calculate expected values for the configuration of the MRT
	// kernel equivalent to LBGK.
	distribn_t expectedDensity0;
	VEC expectedMomentum0;
	DISTS expectedFEq0;
	LbTestsHelper::CalculateRhoMomentum<LATTICE>(hydroVars0.f, expectedDensity0, expectedMomentum0);
	LbTestsHelper::CalculateLBGKEqmF<LATTICE>(expectedDensity0,
						  expectedMomentum0,
						  expectedFEq0);

	// Now compare the expected and actual values.
	LbTestsHelper::CompareHydros(expectedDensity0,
				     expectedMomentum0,
				     expectedFEq0,
				     hydroVars0,
				     allowedError);

	// Do the MRT collision.
	mrtLbgkEquivalentKernel.DoCollide(lbmParams, hydroVars0);

	// Get the expected post-collision velocity distributions with LBGK.
	DISTS expectedPostCollision0;
	LbTestsHelper::CalculateLBGKCollision<LATTICE>(f_original,
						       hydroVars0.GetFEq().f,
						       lbmParams->GetOmega(),
						       expectedPostCollision0);

	// Compare.
	for (unsigned int ii = 0; ii < NV; ++ii) {
	  REQUIRE(Approx(expectedPostCollision0[ii]).margin(allowedError) == hydroVars0.GetFPostCollision()[ii]);
	}
      }
    };
    using fix15 = MRTTestFixture<lb::lattices::D3Q15, lb::kernels::momentBasis::DHumieresD3Q15MRTBasis>;
    TEST_CASE_METHOD(fix15, "MRT with constant relaxation time equals LBGK (15 velocity)") {
    }
    using fix19 = MRTTestFixture<lb::lattices::D3Q19, lb::kernels::momentBasis::DHumieresD3Q19MRTBasis>;
    TEST_CASE_METHOD(fix19, "MRT with constant relaxation time equals LBGK (19 velocity)") {
    }

  }
}

