#ifndef HEMELBSETUPTOOL_TEST_SURFACEVOXELISER_HPP
#define HEMELBSETUPTOOL_TEST_SURFACEVOXELISER_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "TestResources/Meshes.hpp"
#include "TriangleSorter.h"
#include "SurfaceVoxeliser.h"
#include "Neighbours.h"

#include "range.hpp"
#include "enumerate.hpp"

class SurfaceVoxeliserTests: public CppUnit::TestFixture {
	CPPUNIT_TEST_SUITE (SurfaceVoxeliserTests);

	CPPUNIT_TEST (NeighbourInverses);
	CPPUNIT_TEST (Trivial);
	CPPUNIT_TEST (Sphere);

	CPPUNIT_TEST_SUITE_END();

public:
	void NeighbourInverses() {
		auto& vectors = Neighbours::GetDisplacements();
		auto& opps = Neighbours::GetOpposites();
		CPPUNIT_ASSERT(vectors.size() == 26);
		CPPUNIT_ASSERT(opps.size() == 26);
		for (auto i : range(26)) {
			CPPUNIT_ASSERT(vectors[i] + vectors[opps[i]] == Index::Zero());
		}
	}

	void Trivial() {
		// 16 cube
		auto levels = 4;
		auto tri_level = 2;
		auto triv = SimpleMeshFactory::MkTrivial();
		auto tree = TrianglesToTreeSerial(levels, tri_level, triv->points,
				triv->triangles);

		SurfaceVoxeliser voxer(1 << tri_level, triv->points, triv->triangles, triv->normals,
				triv->labels);
		auto tri_node = tree.Get(0,0,0, tri_level);
		auto edge_node = voxer.ComputeIntersectionsForRegion(tri_node);

		auto dirs = Neighbours::GetDisplacements();
		edge_node->IterDepthFirst(0, 0, [&](VoxTree::NodePtr node) {
			// the 2 triangles are at x = 1.2
			// with y = {1.2, 2.2}
			//  and z = {1.2, 2.2}
				CPPUNIT_ASSERT(node->X() == 1 || node->X() == 2);
				auto cuts = node->Data()->closest_cut;
				for (int i = 0; i<26; ++i)
				if (cuts[i].id >= 0) {
					double expected = dirs[i].x > 0 ? 0.2 : 0.8;
					CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, cuts[i].dist, 1e-9);
				}
			});
	}

	void Sphere() {
		TriTree::Int levels = 5;
		TriTree::Int tri_level = 3;

		Vector sphere_centre(10.5);
		double sphere_radius = 10.0;

		auto sphere = SimpleMeshFactory::MkSphere();
		auto tree = TrianglesToTreeSerial(levels, tri_level, sphere->points,
				sphere->triangles);

		SurfaceVoxeliser voxer(1 << tri_level, sphere->points, sphere->triangles,
				sphere->normals, sphere->labels);
		auto edge_tree = voxer(tree, tri_level);
		// Assert things...
		edge_tree.IterDepthFirst([&](VoxTree::NodePtr node) {
			// Has a valid data pointer iff a leaf node
			if (node->Level() == 0) {
				CPPUNIT_ASSERT(node->Data());
			} else {
				CPPUNIT_ASSERT(!node->Data());
				// No further tests for non-leaf
				return;
			}
			Index coord(node->X(), node->Y(), node->Z());
			const auto& cuts = node->Data()->closest_cut;
			int nCuts = 0;
			for (auto& c: cuts)
				nCuts += c.dist < 1.0 ? 1 : 0;

			// All voxels must have 1 or more cuts
			CPPUNIT_ASSERT(nCuts > 0);
		});

		// 24, 17, 20 is an arbitrary exterior point near the surface
		auto vox = edge_tree.Get(24, 17, 20, 0);
		CPPUNIT_ASSERT(vox);
		Index coord(24, 17, 20);

		const auto& cuts = vox->Data()->closest_cut;
		auto directions = Neighbours::GetDisplacements();
		for (auto i: range(26)) {
			const Index& dir = directions[i];
			switch (dir.x) {
			case 1:
				// All +x vectors must have no cut
				CPPUNIT_ASSERT(cuts[i].dist > 1.0);
				break;
			case 0:
				// x = our point, O = outside, I = inside
				// OOO
				// OxO
				// IIO
				if (dir.z == -1 && (dir.y != 1)) {
					// only these two inside
					CPPUNIT_ASSERT(cuts[i].dist < 1.0);
					CPPUNIT_ASSERT(cuts[i].id == 18);
				} else {
					// other 6 outside
					CPPUNIT_ASSERT(cuts[i].dist > 1.0);
				}
				break;
			case -1:
				// x = our point, O = outside, I = inside
				// I*O
				// IxI
				// III
				// * == inside but crossing tri 19
				if (dir.y ==1 && dir.z ==1) {
					// outside
					CPPUNIT_ASSERT(cuts[i].dist > 1.0);
				} else {
					// inside
					CPPUNIT_ASSERT(cuts[i].dist < 1.0);
					if (dir == Index(-1, 0, 1)) {
						CPPUNIT_ASSERT(cuts[i].id == 19);
					} else {
						CPPUNIT_ASSERT(cuts[i].id == 18);
					}
				}
				break;
			default:
				CPPUNIT_FAIL("Invalid direction");
			}
		}
	}
};

CPPUNIT_TEST_SUITE_REGISTRATION (SurfaceVoxeliserTests);

#endif
