// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <iostream>
#include <sstream>

#include <catch2/catch.hpp>

#include "lb/kernels/Kernels.h"
#include "lb/streamers/Streamers.h"
#include "geometry/SiteData.h"

#include "tests/helpers/FourCubeBasedTestFixture.h"
#include "tests/lb/LbTestsHelper.h"

namespace hemelb
{
  namespace tests
  {
    constexpr distribn_t allowedError = 1e-10;

    // StreamerTests
    //
    // This class tests the streamer implementations. We assume the
    // collision operators are correct (as they're tested elsewhere),
    // then compare the post-streamed values with the values we expect
    // to have been streamed there.
    
    TEST_CASE_METHOD(public helpers::FourCubeBasedTestFixture<>, "StreamerTests") {
      using LATTICE = lb::lattices::D3Q15;
      using KERNEL = lb::kernels::LBGK<LATTICE>;
      using COLLISION = lb::collisions::Normal<KERNEL>;
      constexpr auto NUMVECTORS = LATTICE::NUMVECTORS;
      auto propertyCache = std::make_unique<lb::MacroscopicPropertyCache>(*simState, *latDat);
      auto normalCollision = std::make_unique<COLLISION>(initParams);
      auto apprx = [&](double x) {
	return Approx(x).margin(allowedError);
      };

      SECTION("SimpleCollideAndStream") {
	lb::streamers::SimpleCollideAndStream<COLLISION> simpleCollideAndStream(initParams);

	// Initialise fOld in the lattice data. We choose values so
	// that each site has an anisotropic distribution function,
	// and that each site's function is distinguishable.
	LbTestsHelper::InitialiseAnisotropicTestData<LATTICE>(latDat);

	// Use the streaming operator on the entire lattice.
	simpleCollideAndStream.StreamAndCollide<false> (0,
							latDat->GetLocalFluidSiteCount(),
							lbmParams,
							latDat,
							*propertyCache);

	// Now, go over each lattice site and check each value in f_new is correct.
	for (site_t streamedToSite = 0; streamedToSite < latDat->GetLocalFluidSiteCount(); ++streamedToSite) {
	  auto streamedSite = latDat->GetSite(streamedToSite);

	  distribn_t* streamedToFNew = latDat->GetFNew(NUMVECTORS * streamedToSite);

	  for (auto streamedDirection = 0U; streamedDirection < NUMVECTORS; ++streamedDirection) {

	    site_t streamerIndex = streamedSite.GetStreamedIndex<LATTICE>(LATTICE::INVERSEDIRECTIONS[streamedDirection]);

	    // If this site streamed somewhere sensible, it must have been streamed to.
	    if (streamerIndex >= 0 && streamerIndex < (NUMVECTORS * latDat->GetLocalFluidSiteCount())) {
	      site_t streamerSiteId = streamerIndex / NUMVECTORS;

	      // Calculate streamerFOld at this site.
	      distribn_t streamerFOld[NUMVECTORS];
	      LbTestsHelper::InitialiseAnisotropicTestData<LATTICE>(streamerSiteId,
								    streamerFOld);

	      // Calculate what the value streamed to site streamedToSite should be.
	      lb::kernels::HydroVars<KERNEL> streamerHydroVars(streamerFOld);
	      normalCollision->CalculatePreCollision(streamerHydroVars, streamedSite);

	      normalCollision->Collide(lbmParams, streamerHydroVars);

	      // F_new should be equal to the value that was streamed
	      // from this other site in the same direction as we're
	      // streaming from.
	      REQUIRE(apprx(streamerHydroVars.GetFPostCollision()[streamedDirection]) == streamedToFNew[streamedDirection]);
	    }
	  }
	}
      }

      SECTION("BouzidiFirdaousLallemand") {
	// Initialise fOld in the lattice data. We choose values so
	// that each site has an anisotropic distribution function,
	// and that each site's function is distinguishable.
	LbTestsHelper::InitialiseAnisotropicTestData<LATTICE>(latDat);
	lb::streamers::BouzidiFirdaousLallemand<COLLISION>::Type bfl(initParams);

	bfl.StreamAndCollide<false> (0,
				     latDat->GetLocalFluidSiteCount(),
				     lbmParams,
				     latDat,
				     *propertyCache);
	bfl.PostStep<false> (0,
			     latDat->GetLocalFluidSiteCount(),
			     lbmParams,
			     latDat,
			     *propertyCache);

	// Now, go over each lattice site and check each value in f_new is correct.
	for (site_t streamedToSite = 0; streamedToSite < latDat->GetLocalFluidSiteCount(); ++streamedToSite) {
	    const auto streamedSite = latDat->GetSite(streamedToSite);

	    distribn_t* streamedToFNew = latDat->GetFNew(NUMVECTORS * streamedToSite);

	    for (unsigned int streamedDirection = 0; streamedDirection < NUMVECTORS; ++streamedDirection) {
	      unsigned int oppDirection = LATTICE::INVERSEDIRECTIONS[streamedDirection];

	      site_t streamerIndex = streamedSite.GetStreamedIndex<LATTICE>(oppDirection);

	      auto streamerSite = latDat->GetSite(streamerIndex);

	      // If this site streamed somewhere sensible, it must
	      // have been streamed to.
	      if (streamerIndex >= 0 && streamerIndex < (NUMVECTORS * latDat->GetLocalFluidSiteCount())) {
		site_t streamerSiteId = streamerIndex / NUMVECTORS;

		// Calculate streamerFOld at this site.
		distribn_t streamerFOld[NUMVECTORS];
		LbTestsHelper::InitialiseAnisotropicTestData<LATTICE>(streamerSiteId,
								      streamerFOld);

		// Calculate what the value streamed to site streamedToSite should be.
		lb::kernels::HydroVars<KERNEL> streamerHydroVars(streamerFOld);
		normalCollision->CalculatePreCollision(streamerHydroVars, streamerSite);

		normalCollision->Collide(lbmParams, streamerHydroVars);

		// F_new should be equal to the value that was
		// streamed from this other site in the same direction
		// as we're streaming from.
		REQUIRE(apprx(streamerHydroVars.GetFPostCollision()[streamedDirection]) == streamedToFNew[streamedDirection]);
	      } else if (streamedSite.GetSiteType() == geometry::INLET_TYPE
			 || streamedSite.GetSiteType() == geometry::OUTLET_TYPE) {
		// No reason to further test an inlet/outlet site.
		// Pass.
	      } else {
		{
		  INFO("Site: " << streamedToSite
		       << " Direction: " << oppDirection
		       << " Data: " << streamedSite.GetSiteData().GetWallIntersectionData());
		  UNSCOPED_INFO("Expected to find a boundary opposite an unstreamed-to direction.");
		  REQUIRE(streamedSite.HasWall(oppDirection));
		  // Test disabled due to RegressionTests issue, see discussion in #87
		  UNSCOPED_INFO("Expect defined cut distance opposite an unstreamed-to direction ");
		  REQUIRE(streamedSite.GetWallDistance<LATTICE>(oppDirection) != -1.0);
		}
		// To verify the operation of the BFL boundary condition, we'll need:
		// - the distance to the wall * 2
		distribn_t twoQ = 2.0 * streamedSite.GetWallDistance<LATTICE>(oppDirection);

		// - the post-collision distribution at the current site.
		distribn_t streamedToSiteFOld[NUMVECTORS];
		// (initialise it to f_old).
		LbTestsHelper::InitialiseAnisotropicTestData<LATTICE>(streamedToSite,
								      streamedToSiteFOld);

		lb::kernels::HydroVars<KERNEL> hydroVars(streamedToSiteFOld);

		normalCollision->CalculatePreCollision(hydroVars, streamedSite);

		// (find post-collision values using the collision operator).
		normalCollision->Collide(lbmParams, hydroVars);

		// - finally, the post-collision distribution at the
		//   site which is one further away from the wall in
		//   this direction.
		distribn_t awayFromWallFOld[NUMVECTORS];

		site_t awayFromWallIndex = streamedSite.GetStreamedIndex<LATTICE> (streamedDirection) / NUMVECTORS;

		// If there's a valid index in that direction, use BFL
		if (awayFromWallIndex >= 0 && awayFromWallIndex  < latDat->GetLocalFluidSiteCount()) {
		  const auto awayFromWallSite = latDat->GetSite(awayFromWallIndex);

		  // (initialise it to f_old).
		  LbTestsHelper::InitialiseAnisotropicTestData<LATTICE>(awayFromWallIndex,
									awayFromWallFOld);

		  lb::kernels::HydroVars<KERNEL> awayFromWallsHydroVars(awayFromWallFOld);

		  normalCollision->CalculatePreCollision(awayFromWallsHydroVars, awayFromWallSite);

		  // (find post-collision values using the collision operator).
		  normalCollision->Collide(lbmParams, awayFromWallsHydroVars);

		  distribn_t toWallOld = hydroVars.GetFPostCollision()[oppDirection];

		  distribn_t toWallNew = awayFromWallsHydroVars.GetFPostCollision()[oppDirection];

		  distribn_t oppWallOld = hydroVars.GetFPostCollision()[streamedDirection];

		  // The streamed value should be as given below.
		  distribn_t streamed = (twoQ < 1.0)
		    ? (toWallNew + twoQ * (toWallOld - toWallNew))
		    : (oppWallOld + (1. / twoQ) * (toWallOld - oppWallOld));

		  INFO("BouzidiFirdaousLallemand, PostStep: site " << streamedToSite
		       << " direction " << streamedDirection);

		  // Assert that this is the case.
		  REQUIRE(apprx(streamed) == *(latDat->GetFNew(streamedToSite * NUMVECTORS + streamedDirection)));
		} else {
		  // With no valid lattice site, simple bounce-back will be performed.
		  INFO("BouzidiFirdaousLallemand, PostStep by simple bounce-back:"
		       << " site " << streamedToSite
		       << " direction " << streamedDirection);
		  REQUIRE(apprx(hydroVars.GetFPostCollision()[oppDirection]) == *(latDat->GetFNew(streamedToSite * NUMVECTORS + streamedDirection)));
		}
	      }
	    }
	}
      }

      SECTION("TestSimpleBounceBack") {
	// Initialise fOld in the lattice data. We choose values so
	// that each site has an anisotropic distribution function,
	// and that each site's function is distinguishable.
	LbTestsHelper::InitialiseAnisotropicTestData<LATTICE>(latDat);

	site_t firstWallSite = latDat->GetMidDomainCollisionCount(0);
	site_t wallSitesCount = latDat->GetMidDomainCollisionCount(1);

	// Check that the lattice has the expected number of sites
	// labeled as pure wall (otherwise this test is void)
	REQUIRE(site_t(24) == latDat->GetMidDomainCollisionCount(1));

	site_t offset = 0;

	// Mid-Fluid sites use simple collide and stream
	lb::streamers::SimpleCollideAndStream<COLLISION > simpleCollideAndStream(initParams);

	simpleCollideAndStream.StreamAndCollide<false> (offset,
							latDat->GetMidDomainCollisionCount(0),
							lbmParams,
							latDat,
							*propertyCache);
	offset += latDat->GetMidDomainCollisionCount(0);

	// Wall sites use simple bounce back
	lb::streamers::SimpleBounceBack<COLLISION>::Type simpleBounceBack(initParams);

	simpleBounceBack.StreamAndCollide<false> (offset,
						  latDat->GetMidDomainCollisionCount(1),
						  lbmParams,
						  latDat,
						  *propertyCache);
	offset += latDat->GetMidDomainCollisionCount(1);

	// Consider inlet/outlets and their walls as mid-fluid sites
	simpleCollideAndStream.StreamAndCollide<false> (offset,
							latDat->GetLocalFluidSiteCount()
							- offset,
							lbmParams,
							latDat,
							*propertyCache);
	offset += latDat->GetLocalFluidSiteCount() - offset;

	// Sanity check
	REQUIRE(offset == latDat->GetLocalFluidSiteCount());

	// Loop over the wall sites and check whether they got
	// properly streamed on or bounced back depending on where
	// they sit relative to the wall. We ignore mid-Fluid sites
	// since StreamAndCollide was tested before.
	for (site_t wallSiteLocalIndex = 0; wallSiteLocalIndex < wallSitesCount; wallSiteLocalIndex++) {
	  site_t streamedToSite = firstWallSite + wallSiteLocalIndex;
	  const auto streamedSite = latDat->GetSite(streamedToSite);
	  distribn_t* streamedToFNew = latDat->GetFNew(NUMVECTORS * streamedToSite);

	  for (unsigned int streamedDirection = 0; streamedDirection
		 < NUMVECTORS; ++streamedDirection) {
	    unsigned oppDirection = LATTICE::INVERSEDIRECTIONS[streamedDirection];

	    // Index of the site streaming to streamedToSite via direction streamedDirection
	    site_t streamerIndex = streamedSite.GetStreamedIndex<LATTICE> (oppDirection);

	    // Is streamerIndex a valid index?
	    if (streamerIndex >= 0 && streamerIndex < (NUMVECTORS * latDat->GetLocalFluidSiteCount())) {
	      // The streamer index is a valid index in the domain,
	      // therefore stream and collide has happened
	      site_t streamerSiteId = streamerIndex / NUMVECTORS;

	      // Calculate streamerFOld at this site.
	      distribn_t streamerFOld[NUMVECTORS];
	      LbTestsHelper::InitialiseAnisotropicTestData<LATTICE>(streamerSiteId,
								    streamerFOld);

	      // Calculate what the value streamed to site streamedToSite should be.
	      lb::kernels::HydroVars<KERNEL> streamerHydroVars(streamerFOld);
	      normalCollision->CalculatePreCollision(streamerHydroVars, streamedSite);

	      normalCollision->Collide(lbmParams, streamerHydroVars);

	      // F_new should be equal to the value that was streamed
	      // from this other site
	      // in the same direction as we're streaming from.
	      INFO("SimpleCollideAndStream, StreamAndCollide");
	      REQUIRE(apprx(streamerHydroVars.GetFPostCollision()[streamedDirection]) == streamedToFNew[streamedDirection]);
	    } else {
	      // The streamer index shows that no one has streamed to
	      // streamedToSite direction streamedDirection, therefore
	      // bounce back has happened in that site for that
	      // direction

	      // Initialise streamedToSiteFOld with the original data
	      distribn_t streamerToSiteFOld[NUMVECTORS];
	      LbTestsHelper::InitialiseAnisotropicTestData<LATTICE>(streamedToSite,
								    streamerToSiteFOld);
	      lb::kernels::HydroVars<KERNEL> hydroVars(streamerToSiteFOld);
	      normalCollision->CalculatePreCollision(hydroVars, streamedSite);

	      // Simulate post-collision using the collision operator.
	      normalCollision->Collide(lbmParams, hydroVars);

	      // After streaming FNew in a given direction must be f
	      // post-collision in the opposite direction following
	      // collision
	      INFO("Simple bounce-back: site " << streamedToSite << " direction " << streamedDirection);
	      REQUIRE(streamedToFNew[streamedDirection] == apprx(hydroVars.GetFPostCollision()[oppDirection]));
	    }
	  }
	}
      }

      SECTION("GuoZhengShi") {
	lb::streamers::GuoZhengShi<COLLISION>::Type guoZhengShi(initParams);

	for (double assignedWallDistance = 0.4;
	     assignedWallDistance < 1.0;
	     assignedWallDistance += 0.5) {
	  // Initialise fOld in the lattice data. We choose values so
	  // that each site has an anisotropic distribution function,
	  // and that each site's function is distinguishable.
	  LbTestsHelper::InitialiseAnisotropicTestData<LATTICE>(latDat);

	  // Make some fairly arbitrary choices early on.
	  const site_t chosenSite = 0;
	  const auto& streamer = latDat->GetSite(chosenSite);

	  const Direction chosenUnstreamedDirection = 5;
	  const Direction chosenWallDirection = LATTICE::INVERSEDIRECTIONS[chosenUnstreamedDirection];
	  const Direction chosenDoubleWallDirection1 = 7;
	  const Direction chosenDoubleWallDirection2 = 8;

	  // Enforce that there's a boundary in the wall direction.
	  latDat->SetHasWall(chosenSite, chosenWallDirection);
	  latDat->SetHasWall(chosenSite, chosenDoubleWallDirection1);
	  latDat->SetHasWall(chosenSite, chosenDoubleWallDirection2);
	  latDat->SetBoundaryDistance(chosenSite, chosenWallDirection, assignedWallDistance);
	  latDat->SetBoundaryDistance(chosenSite,
				      chosenDoubleWallDirection1,
				      assignedWallDistance);
	  latDat->SetBoundaryDistance(chosenSite,
				      chosenDoubleWallDirection2,
				      assignedWallDistance);

	  // Perform the collision and streaming.
	  guoZhengShi.StreamAndCollide<false> (chosenSite, 1, lbmParams, latDat, *propertyCache);

	  // Calculate the distributions at the chosen site up to post-collision.
	  distribn_t streamerFOld[NUMVECTORS];
	  LbTestsHelper::InitialiseAnisotropicTestData<LATTICE>(chosenSite,
								streamerFOld);
	  lb::kernels::HydroVars<KERNEL> streamerHydroVars(streamerFOld);
	  normalCollision->CalculatePreCollision(streamerHydroVars, streamer);
	  normalCollision->Collide(lbmParams, streamerHydroVars);

	  // Check each streamed direction.
	  for (Direction streamedDirection = 0; streamedDirection < NUMVECTORS; ++streamedDirection) {
	    switch (streamedDirection) {
	    case 5:
	      {
		// We point away from the wall and there is a fluid
		// site out in this direction, so we have regular GZS.

		LatticeVelocity velocityWall;
		distribn_t fNeqWall;

		// This is the first means of estimating from the
		// source paper: only use the nearest fluid site.
		LatticeVelocity velocityEstimate1 = streamerHydroVars.momentum *
		  ((1. - 1./ assignedWallDistance) / streamerHydroVars.density);

		distribn_t fNeqEstimate1 = streamerHydroVars.GetFNeq()[streamedDirection];

		if (assignedWallDistance < 0.75) {
		  // This is the second method for estimating: using
		  // the next fluid site away from the wall.
		  const site_t nextSiteAwayFromWall = streamer.GetStreamedIndex<LATTICE> (streamedDirection) / NUMVECTORS;
		  const auto& nextSiteAway = latDat->GetSite(nextSiteAwayFromWall);
		  distribn_t nextSiteOutFOld[NUMVECTORS];
		  LbTestsHelper::InitialiseAnisotropicTestData<LATTICE>(nextSiteAwayFromWall,
									nextSiteOutFOld);
		  lb::kernels::HydroVars<KERNEL> nextSiteOutHydroVars(nextSiteOutFOld);
		  normalCollision->CalculatePreCollision(nextSiteOutHydroVars, nextSiteAway);

		  LatticeVelocity velocityEstimate2 = nextSiteOutHydroVars.momentum
		    * ((assignedWallDistance - 1.) /
		       ( (1. + assignedWallDistance) * nextSiteOutHydroVars.density));

		  distribn_t fNeqEstimate2 = nextSiteOutHydroVars.GetFNeq()[streamedDirection];

		  // The actual value is taken to be an interpolation
		  // between the two estimates.
		  velocityWall = velocityEstimate1 * assignedWallDistance
		    + velocityEstimate2 * (1. - assignedWallDistance);

		  fNeqWall = assignedWallDistance * fNeqEstimate1
		    + (1. - assignedWallDistance) * fNeqEstimate2;
		} else {
		  // Only use the first estimate
		  velocityWall = velocityEstimate1;
		  fNeqWall = fNeqEstimate1;
		}

		auto momentumWall = velocityWall * streamerHydroVars.density;

		// Compute the wall equilibrium dist
		distribn_t fEqm[NUMVECTORS];
		LATTICE::CalculateFeq(streamerHydroVars.density,
				      momentumWall.x,
				      momentumWall.y,
				      momentumWall.z,
				      fEqm);
		// Perform collision on the wall f's
		distribn_t prediction = fEqm[streamedDirection] + (1.0 + lbmParams->GetOmega()) * fNeqWall;
		// This is the answer from the code we're testing
		distribn_t streamedFNew = latDat->GetFNew(NUMVECTORS * chosenSite)[streamedDirection];
		REQUIRE(apprx(prediction) == streamedFNew);
		break;
	      }
	    case 7:
	    case 8:
	      // These direction both point at the wall and are
	      // opposite each other.  We might do SBB or GZS
	      // depending on how close the wall is.
	      if (assignedWallDistance < 0.75) {
		// It's SBB
		Direction inv = LATTICE::INVERSEDIRECTIONS[streamedDirection];
		distribn_t prediction = streamerHydroVars.GetFPostCollision()[inv];
		// This is the answer from the code we're testing
		distribn_t streamedFNew = latDat->GetFNew(NUMVECTORS * chosenSite)[streamedDirection];
		REQUIRE(apprx(prediction) == streamedFNew);
	      } else {
		// It's GZS with extrapolation from this site only

		// This is the first means of estimating from the
		// source paper: only use the nearest fluid site.
		LatticeVelocity velocityWall = streamerHydroVars.momentum * ((1. - 1./ assignedWallDistance) / streamerHydroVars.density);

		distribn_t fNeqWall = streamerHydroVars.GetFNeq()[streamedDirection];

		auto momentumWall = velocityWall * streamerHydroVars.density;

		// Compute the wall equilibrium dist
		distribn_t fEqm[NUMVECTORS];
		LATTICE::CalculateFeq(streamerHydroVars.density,
				      momentumWall.x,
				      momentumWall.y,
				      momentumWall.z,
				      fEqm);

		// Perform collision on the wall f's
		distribn_t prediction = fEqm[streamedDirection] + (1.0 + lbmParams->GetOmega()) * fNeqWall;
		// This is the answer from the code we're testing
		distribn_t streamedFNew = latDat->GetFNew(NUMVECTORS * chosenSite)[streamedDirection];
		REQUIRE(apprx(prediction) == streamedFNew);
	      }
	      break;
	    case 6:
	      // We point towards the wall.
	      break;

	    default:
	      // We have nothing to do with a wall so simple streaming
	      const site_t streamedIndex = streamer.GetStreamedIndex<LATTICE> (streamedDirection);
	      distribn_t streamedToFNew = *latDat->GetFNew(streamedIndex);

	      // F_new should be equal to the value that was streamed
	      // from this other site in the same direction as we're
	      // streaming from.
	      REQUIRE(apprx(streamerHydroVars.GetFPostCollision()[streamedDirection]) == streamedToFNew);
	      break;
	    }
	  }
	}
      }

      // Junk&Yang should behave like simple bounce back when fluid
      // sites are 0.5 lattice length units away from the domain
      // boundary.
      SECTION("JunkYangEquivalentToBounceBack") {
	// Initialise fOld in the lattice data. We choose values so
	// that each site has an anisotropic distribution function,
	// and that each site's function is distinguishable.
	LbTestsHelper::InitialiseAnisotropicTestData<LATTICE>(latDat);
	// Setting all the wall distances to 0.5 will make Junk&Yang
	// behave like Simple Bounce Back
	LbTestsHelper::SetWallAndIoletDistances<LATTICE>(*latDat, 0.5);

	site_t firstWallSite = latDat->GetMidDomainCollisionCount(0);
	site_t wallSitesCount = latDat->GetMidDomainCollisionCount(1) - firstWallSite;

	// Check that the lattice has sites labeled as wall (otherwise
	// this test is void)
	REQUIRE(wallSitesCount == 16);

	site_t offset = 0;

	// Mid-Fluid sites use simple collide and stream
	lb::streamers::SimpleCollideAndStream<COLLISION> simpleCollideAndStream(initParams);

	simpleCollideAndStream.StreamAndCollide<false> (offset,
							latDat->GetMidDomainCollisionCount(0),
							lbmParams,
							latDat,
							*propertyCache);
	offset += latDat->GetMidDomainCollisionCount(0);

	// Wall sites use Junk and Yang
	initParams.siteRanges.push_back(std::pair<site_t, site_t>(offset,
								  offset
								  + latDat->GetMidDomainCollisionCount(1)));
	lb::streamers::JunkYang<COLLISION>::Type junkYang(initParams);

	junkYang.StreamAndCollide<false> (offset,
					  latDat->GetMidDomainCollisionCount(1),
					  lbmParams,
					  latDat,
					  *propertyCache);

	junkYang.PostStep<false> (offset,
				  latDat->GetMidDomainCollisionCount(1),
				  lbmParams,
				  latDat,
				  *propertyCache);

	offset += latDat->GetMidDomainCollisionCount(1);

	// Consider inlet/outlets and their walls as mid-fluid sites
	simpleCollideAndStream.StreamAndCollide<false> (offset,
							latDat->GetLocalFluidSiteCount()
							- offset,
							lbmParams,
							latDat,
							*propertyCache);
	offset += latDat->GetLocalFluidSiteCount() - offset;

	// Sanity check
	REQUIRE(offset == latDat->GetLocalFluidSiteCount());

	// Loop over the wall sites and check whether they got
	// properly streamed on or bounced back depending on where
	// they sit relative to the wall. We ignore mid-Fluid sites
	// since StreamAndCollide was tested before.
	for (site_t wallSiteLocalIndex = 0; wallSiteLocalIndex < wallSitesCount; wallSiteLocalIndex++) {
	  site_t streamedToSite = firstWallSite + wallSiteLocalIndex;
	  const auto streamedSite = latDat->GetSite(streamedToSite);
	  distribn_t* streamedToFNew = latDat->GetFNew(NUMVECTORS * streamedToSite);

	  for (unsigned int streamedDirection = 0;
	       streamedDirection < NUMVECTORS; ++streamedDirection) {
	    unsigned oppDirection = LATTICE::INVERSEDIRECTIONS[streamedDirection];

	    // Index of the site streaming to streamedToSite via
	    // direction streamedDirection
	    site_t streamerIndex = streamedSite.GetStreamedIndex<LATTICE> (oppDirection);

	    // Is streamerIndex a valid index?
	    if (streamerIndex >= 0 &&
		streamerIndex < (NUMVECTORS * latDat->GetLocalFluidSiteCount())) {
	      // The streamer index is a valid index in the domain,
	      // therefore stream and collide has happened
	      site_t streamerSiteId = streamerIndex / NUMVECTORS;

	      // Calculate streamerFOld at this site.
	      distribn_t streamerFOld[NUMVECTORS];
	      LbTestsHelper::InitialiseAnisotropicTestData<LATTICE>(streamerSiteId,
								    streamerFOld);

	      // Calculate what the value streamed to site streamedToSite should be.
	      lb::kernels::HydroVars<KERNEL> streamerHydroVars(streamerFOld);
	      normalCollision->CalculatePreCollision(streamerHydroVars, streamedSite);

	      normalCollision->Collide(lbmParams, streamerHydroVars);

	      // F_new should be equal to the value that was streamed
	      // from this other site in the same direction as we're
	      // streaming from.
	      REQUIRE(apprx(streamerHydroVars.GetFPostCollision()[streamedDirection])
		      == streamedToFNew[streamedDirection]);
	    } else {
	      // The streamer index shows that no one has streamed to
	      // streamedToSite direction streamedDirection, therefore
	      // bounce back has happened in that site for that
	      // direction

	      // Initialise streamedToSiteFOld with the original data
	      distribn_t streamerToSiteFOld[NUMVECTORS];
	      LbTestsHelper::InitialiseAnisotropicTestData<LATTICE>(streamedToSite, streamerToSiteFOld);
	      lb::kernels::HydroVars<KERNEL> hydroVars(streamerToSiteFOld);
	      normalCollision->CalculatePreCollision(hydroVars, streamedSite);

	      // Simulate post-collision using the collision operator.
	      normalCollision->Collide(lbmParams, hydroVars);

	      // After streaming FNew in a given direction must be f
	      // post-collision in the opposite direction following
	      // collision
	      INFO("Junk&Yang bounce-back equivalent: site " << streamedToSite
		   << " direction " << streamedDirection);
	      REQUIRE(apprx(streamedToFNew[streamedDirection])
		      == hydroVars.GetFPostCollision()[oppDirection]);
	    }
	  }
	}
      }

      SECTION("NashZerothOrderPressureIolet") {
	lb::iolets::BoundaryValues inletBoundary(geometry::INLET_TYPE,
						 latDat,
						 simConfig->GetInlets(),
						 simState.get(),
						 Comms(),
						 *unitConverter);

	initParams.boundaryObject = &inletBoundary;

	lb::streamers::NashZerothOrderPressureIolet<COLLISION>::Type ioletCollider(initParams);
	    
	for (double assignedWallDistance = 0.4;
	     assignedWallDistance < 1.0;
	     assignedWallDistance += 0.5) {
	  // Initialise fOld in the lattice data. We choose values so
	  // that each site has an anisotropic distribution function,
	  // and that each site's function is distinguishable.
	  LbTestsHelper::InitialiseAnisotropicTestData<LATTICE>(latDat);

	  // Make some fairly arbitrary choices early on.
	  const site_t chosenSite = 0;
	  const int chosenBoundaryId = 0;
	  const auto& streamer = latDat->GetSite(chosenSite);

	  const Direction chosenUnstreamedDirection = 5;
	  const Direction chosenIoletDirection = LATTICE::INVERSEDIRECTIONS[chosenUnstreamedDirection];
	  const auto ioletNormal = inletBoundary.GetLocalIolet(chosenBoundaryId)->GetNormal();

	  // Enforce that there's a boundary in the iolet direction.
	  latDat->SetHasIolet(chosenSite, chosenIoletDirection);
	  latDat->SetBoundaryDistance(chosenSite, chosenIoletDirection, assignedWallDistance);
	  latDat->SetBoundaryNormal(chosenSite, ioletNormal);
	  latDat->SetIoletId(chosenSite, chosenBoundaryId);

	  // Perform the collision and streaming.
	  ioletCollider.StreamAndCollide<false> (chosenSite,
						 1,
						 lbmParams,
						 latDat,
						 *propertyCache);

	  // Check each streamed direction.
	  for (Direction streamedDirection = 0; streamedDirection < NUMVECTORS; ++streamedDirection) {
	    // Calculate the distributions at the chosen site up to post-collision.
	    distribn_t streamerFOld[NUMVECTORS];
	    LbTestsHelper::InitialiseAnisotropicTestData<LATTICE>(chosenSite,
								  streamerFOld);

	    lb::kernels::HydroVars<KERNEL> streamerHydroVars(streamerFOld);
	    normalCollision->CalculatePreCollision(streamerHydroVars, streamer);
	    normalCollision->Collide(lbmParams, streamerHydroVars);

	    // Calculate the streamed-to index.
	    const site_t streamedIndex = streamer.GetStreamedIndex<LATTICE> (streamedDirection);

	    // Check that simple collide and stream has happened when
	    // appropriate.  Is streamerIndex a valid index? (And is
	    // it not in one of the directions that has been meddled
	    // with for the test)?
	    if (!streamer.HasIolet(streamedDirection)
		&& streamedIndex >= 0
		&& streamedIndex < (NUMVECTORS * latDat->GetLocalFluidSiteCount())) {
	      distribn_t streamedToFNew = *latDat->GetFNew(streamedIndex);

	      // F_new should be equal to the value that was streamed
	      // from this other site in the same direction as we're
	      // streaming from.
	      REQUIRE(streamerHydroVars.GetFPostCollision()[streamedDirection]
		      == apprx(streamedToFNew));
	    }

	    // Next, handle the case where this is the direction where
	    // we're checking for behaviour with a wall. I.e. are we
	    // correctly filling distributions that aren't streamed-to
	    // by simple streaming?
	    if (streamedDirection == chosenUnstreamedDirection) {
	      // The streamer works by assuming the presence of a
	      // 'ghost' site, just beyond the iolet. The density of
	      // the ghost site is extrapolated from the iolet density
	      // and the density of the fluid site.
	      distribn_t ghostSiteDensity = inletBoundary.GetBoundaryDensity(chosenBoundaryId);
		
	      // The velocity of the ghost site is the component of
	      // the fluid site's velocity along the iolet normal.
	      auto ghostSiteVelocity = ioletNormal * (streamerHydroVars.momentum / streamerHydroVars.density).Dot(ioletNormal);

	      auto ghostSiteMomentum = ghostSiteVelocity * ghostSiteDensity;

	      distribn_t ghostPostCollision[NUMVECTORS];

	      LbTestsHelper::CalculateLBGKEqmF<LATTICE>(ghostSiteDensity,
							ghostSiteMomentum,
							ghostPostCollision);

	      REQUIRE(latDat->GetFNew(chosenSite * NUMVECTORS)[chosenUnstreamedDirection]
		      == apprx(ghostPostCollision[chosenUnstreamedDirection]));
	    }
	  }
	}
      }

      SECTION("NashZerothOrderPressureBB") {
	lb::iolets::BoundaryValues inletBoundary(geometry::INLET_TYPE,
						 latDat,
						 simConfig->GetInlets(),
						 simState.get(),
						 Comms(),
						 *unitConverter);

	initParams.boundaryObject = &inletBoundary;

	lb::streamers::NashZerothOrderPressureIoletSBB<COLLISION>::Type ioletCollider(initParams);

	for (double assignedIoletDistance = 0.4;
	     assignedIoletDistance < 1.0;
	     assignedIoletDistance += 0.5) {
	  // Initialise fOld in the lattice data. We choose values so
	  // that each site has an anisotropic distribution function,
	  // and that each site's function is distinguishable.
	  LbTestsHelper::InitialiseAnisotropicTestData<LATTICE>(latDat);

	  // Make some fairly arbitrary choices early on.
	  const site_t chosenSite = 0;
	  const int chosenBoundaryId = 0;
	  const geometry::Site<geometry::LatticeData>& streamer = latDat->GetSite(chosenSite);

	  const Direction chosenWallDirection = 11;
	  const Direction chosenUnstreamedDirection = 5;
	  const Direction chosenIoletDirection = LATTICE::INVERSEDIRECTIONS[chosenUnstreamedDirection];
	  const auto ioletNormal = inletBoundary.GetLocalIolet(chosenBoundaryId)->GetNormal();

	  // Enforce that there's a boundary in the iolet direction.
	  latDat->SetHasIolet(chosenSite, chosenIoletDirection);
	  latDat->SetHasWall(chosenSite, chosenWallDirection);
	  latDat->SetBoundaryDistance(chosenSite, chosenIoletDirection, assignedIoletDistance);
	  latDat->SetBoundaryNormal(chosenSite, ioletNormal);
	  latDat->SetIoletId(chosenSite, chosenBoundaryId);

	  // Perform the collision and streaming.
	  ioletCollider.StreamAndCollide<false> (chosenSite,
						 1,
						 lbmParams,
						 latDat,
						 *propertyCache);

	  // Check each streamed direction.
	  for (Direction streamedDirection = 0;
	       streamedDirection < NUMVECTORS; ++streamedDirection) {
	    // Calculate the distributions at the chosen site up to post-collision.
	    distribn_t streamerFOld[NUMVECTORS];
	    LbTestsHelper::InitialiseAnisotropicTestData<LATTICE>(chosenSite,
								  streamerFOld);
	    
	    lb::kernels::HydroVars<KERNEL> streamerHydroVars(streamerFOld);
	    normalCollision->CalculatePreCollision(streamerHydroVars, streamer);
	    normalCollision->Collide(lbmParams, streamerHydroVars);

	    // Calculate the streamed-to index.
	    const site_t streamedIndex = streamer.GetStreamedIndex<LATTICE> (streamedDirection);

	    // Check that simple collide and stream has happened when
	    // appropriate.  Is streamerIndex a valid index? (And is
	    // it not in one of the directions that has been meddled
	    // with for the test)?
	    if (!streamer.HasIolet(streamedDirection)
		&& !streamer.HasWall(streamedDirection)
		&& streamedIndex >= 0
		&& streamedIndex < (NUMVECTORS * latDat->GetLocalFluidSiteCount())) {
	      distribn_t streamedToFNew = *latDat->GetFNew(streamedIndex);

	      // F_new should be equal to the value that was streamed
	      // from this other site in the same direction as we're
	      // streaming from.
	      REQUIRE(apprx(streamerHydroVars.GetFPostCollision()[streamedDirection])
		      == streamedToFNew);
	    }

	    Direction inverseDirection = LATTICE::INVERSEDIRECTIONS[streamedDirection];

	    // Check the case by a wall.
	    if (streamer.HasWall(streamedDirection)) {
	      distribn_t streamedToFNew = *latDat->GetFNew(NUMVECTORS * chosenSite + inverseDirection);

	      REQUIRE(apprx(streamerHydroVars.GetFPostCollision()[streamedDirection])
		      == streamedToFNew);
	    }

	    // Next, handle the case where this is the direction where
	    // we're checking for behaviour with a iolet. I.e. are we
	    // correctly filling distributions that aren't streamed-to
	    // by simple streaming?
	    if (streamedDirection == chosenUnstreamedDirection) {
	      // The streamer works by assuming the presence of a
	      // 'ghost' site, just beyond the iolet. The density of
	      // the ghost site is extrapolated from the iolet density
	      // and the density of the fluid site.
	      distribn_t ghostSiteDensity = inletBoundary.GetBoundaryDensity(chosenBoundaryId);

	      // The velocity of the ghost site is the component of
	      // the fluid site's velocity along the iolet normal.
	      auto ghostSiteVelocity = ioletNormal * (streamerHydroVars.momentum / streamerHydroVars.density).Dot(ioletNormal);

	      auto ghostSiteMomentum = ghostSiteVelocity * ghostSiteDensity;

	      distribn_t ghostPostCollision[NUMVECTORS];

	      LbTestsHelper::CalculateLBGKEqmF<LATTICE>(ghostSiteDensity,
							ghostSiteMomentum,
							ghostPostCollision);

	      REQUIRE(latDat->GetFNew(chosenSite * NUMVECTORS)[chosenUnstreamedDirection]
		      == apprx(ghostPostCollision[chosenUnstreamedDirection]));
	    }
	  }
	}
      }
    }
  }
}
