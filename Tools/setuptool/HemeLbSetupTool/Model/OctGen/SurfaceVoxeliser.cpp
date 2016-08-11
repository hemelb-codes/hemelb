#include <cassert>
#include "SurfaceVoxeliser.h"
#include "MkCgalMesh.h"
#include "Neighbours.h"
#include "SegmentFactory.h"

SurfaceVoxeliser::SurfaceVoxeliser(const int ns,
		const std::vector<Vector>& p,
		const std::vector<Index>& t, const std::vector<Vector>&n,
		const std::vector<int>& l) :
		Points(p), Triangles(t), Normals(n), Labels(l),
		mesh(MkCgalMesh(p, t, l)),
		searcher(new CgalSearchTree(faces(*mesh).first, faces(*mesh).second, *mesh))
{
	SegmentFactory cube(ns);
	const auto& opps = Neighbours::GetOpposites();
	const auto& neighs = Neighbours::GetDisplacements();

	direction_indices.reserve(NDIR);
	directions.reserve(NDIR);
	relative_segments.reserve(NDIR);

	for (int i = 0; i < 2*NDIR; ++i) {
		int iOpp = opps[i];
		if (iOpp > i) {
			direction_indices.push_back(i);
			directions.push_back(neighs[i]);
			relative_segments.emplace_back(cube.MakeSegments(neighs[i]));
		}
	}
	assert(relative_segments.size() == NDIR);
}

void AddIntersection(VoxTree::NodePtr ans,
		const Index& idx, int iDir, double cut_dist, int id) {
	try {
		auto vox = ans->GetCreate(idx.x, idx.y, idx.z, 0);
		if (!vox->Data())
			vox->Data() = std::make_shared<Vox>();
		auto& cuts = vox->Data()->closest_cut;
		if (cut_dist < cuts[iDir].dist) {
			cuts[iDir].dist = cut_dist;
			cuts[iDir].id = id;
		}
	} catch (std::out_of_range& e) {
		// this is on a neighbouring subtree
	}
}

VoxTree::NodePtr SurfaceVoxeliser::ComputeIntersectionsForRegion(const TriTree::Node& node) {
	// We're going to work on lines that span the whole node,
	// starting from and ending at the halo voxels.
	typedef boost::optional< CgalSearchTree::Intersection_and_primitive_id<CgalSegment>::Type > Segment_intersection;

	VoxTree::NodePtr ans(new VoxTree::Branch(node.X(), node.Y(), node.Z(), node.Level()));
	Index node_origin(node.X(), node.Y(), node.Z());

	for (int i = 0; i < NDIR; ++i) {
		Index direction = directions[i];
		int iDir = direction_indices[i];
		int iOpp = Neighbours::GetOpposites()[iDir];

		double norm = std::sqrt(direction.GetMagnitudeSquared());
		for (const auto& seg: relative_segments[i]) {
			auto startIdx = node_origin + seg.first;
			auto endIdx = node_origin + seg.second;
			CgalPoint start(startIdx.x, startIdx.y, startIdx.z);
			CgalPoint end(endIdx.x, endIdx.y, endIdx.z);
			CgalSegment link(start, end);
			// compute all intersections with segment query (as pairs object - primitive_id)
			std::list<Segment_intersection> intersections;
			searcher->all_intersections(link, std::back_inserter(intersections));

			for (auto isec: intersections) {
				CgalPoint& pt = boost::get<CgalPoint>(isec->first);
				auto id = isec->second->id();
				auto link_dist = std::sqrt(CGAL::squared_distance(start, pt)) / norm;
				int offset = link_dist;
				double cut_dist = link_dist - offset;

				auto idxA = startIdx + direction*offset;
				auto idxB = startIdx + direction*(offset+1);
				AddIntersection(ans, idxA, iDir, cut_dist, id);
				AddIntersection(ans, idxB, iOpp, 1-cut_dist, id);
			}

		}
	}
	return ans;
}
