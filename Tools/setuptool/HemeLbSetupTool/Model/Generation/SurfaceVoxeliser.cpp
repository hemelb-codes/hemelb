// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#include <cassert>
#include <deque>

#include "SurfaceVoxeliser.h"
#include "MkCgalMesh.h"
#include "Neighbours.h"
#include "SegmentFactory.h"
#include "ParallelApply.h"
#include "range.hpp"

SurfaceVoxeliser::SurfaceVoxeliser(const int ns,
				   const std::vector<Vector>& p,
				   const std::vector<Index>& t, const std::vector<Vector>&n,
				   const std::vector<int>& l, const std::vector<Iolet>&io) :
  Points(p), Triangles(t), Normals(n), Labels(l), Iolets(io),
  mesh(MkCgalMesh(p, t)),
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

void AddIntersection(EdgeSiteTree::NodePtr ans,
		     const Index& idx, int iDir, double cut_dist, int id) {
  try {
    auto vox = ans->GetCreate(idx.x, idx.y, idx.z, 0);
    if (!vox->Data())
      vox->Data() = std::make_shared<EdgeSite>();
    auto& cuts = vox->Data()->closest_cut;
    if (cut_dist < cuts[iDir].dist) {
      cuts[iDir].dist = cut_dist;
      cuts[iDir].id = id;
    }
  } catch (std::out_of_range& e) {
    // this is on a neighbouring subtree
  }
}

typedef boost::optional< CgalSearchTree::Intersection_and_primitive_id<CgalSegment>::Type > Segment_intersection;
typedef std::pair<double, Segment_intersection> dist_isec;

// Helper class for filtering out clusters of intersections
struct Clusterer {
  const std::vector<Vector>& normals;
  std::list<dist_isec>& output;
  const Vector& direction;
  
  int n_pts;
  int n_pos;
  int n_neg;
  const dist_isec* current;
  
  Clusterer(const std::vector<Vector>& n, std::list<dist_isec>& o, const Vector& d)
    : normals(n), output(o), direction(d), n_pts(0), n_pos(0), n_neg(0), current(nullptr) {
  }

  const double& GetEnd() {
    return current->first;
  }
  
  void Reset() {
    n_pts = 0;
    n_pos = 0;
    n_neg = 0;
    current = nullptr;
  }

  void Start(const dist_isec& di) {
    Reset();
    Extend(di);
    current = &di;
  }
  
  void Extend(const dist_isec& di) {
    n_pts++;
    auto isec = di.second;
    auto face_id = isec->second->id();
    auto dp = Vector::Dot(normals[face_id], direction);
    if (dp > 0.0) {
      n_pos++;
    } else {
      n_neg++;
    }
    current = &di;
  }
  
  void Finish() {
    if (n_pos && n_neg) {
      // do not add anything if 
    } else {
      output.push_back(*current);
    }
  }
};

// This filters a list of intersections for clusters of intersections that should be merged
std::list<dist_isec> FilterClusters(const CgalPoint& start, const Vector direction, const std::vector<Vector>& normals,
				    const std::list<Segment_intersection>& intersections, double tol=1e-3) {
  std::list<dist_isec> ans;
  auto link_dist = [&](const Segment_intersection& isec) -> double {
    const CgalPoint& pt = boost::get<CgalPoint>(isec->first);
    return std::sqrt(CGAL::squared_distance(start, pt));
  };
  
  if (intersections.size() > 1) {
    // iterator pointing to elememt, dist squared, in cluster flag
    std::vector<dist_isec> dists;
    dists.reserve(intersections.size());
    for (auto in_it = intersections.begin(); in_it != intersections.end(); ++in_it) {
      dists.emplace_back(link_dist(*in_it), *in_it);
    }
    
    // Sort by dist squared
    std::sort(dists.begin(), dists.end(),
	      [](const dist_isec& a, const dist_isec& b) {
		return a.first < b.first;
	      });
    
    // Look for "clusters" of points where they are really close together
    
    // Cluster intersections should either be ignored (if a mix of
    // inner and outer) or merged to a single intersection (with an
    // arbitrary face picked).    
    bool first = true;
    Clusterer clst(normals, ans, direction);
    for (const auto& cur: dists) {
      if (first) {
	clst.Start(cur);
	first = false;
	continue;
      }
      
      if (cur.first < (clst.GetEnd() + tol)) {
	clst.Extend(cur);
      } else {
	clst.Finish();
	clst.Start(cur);
      }
    }
    clst.Finish();
    
  } else {
    // zero or one intersections: just compute distances and add to output
    std::for_each(intersections.begin(), intersections.end(),
		  [&ans, &link_dist](const Segment_intersection& isec) {
		    auto d = link_dist(isec);
		    ans.emplace_back(d, isec);
		  });
  }
  return ans;
}

EdgeSiteTree::NodePtr SurfaceVoxeliser::ComputeIntersectionsForRegion(TriTree::ConstNodePtr node) const {
  // We're going to work on lines that span the whole node,
  // starting from and ending at the halo voxels.
  EdgeSiteTree::NodePtr ans(new EdgeSiteTree::Branch(node->X(), node->Y(), node->Z(), node->Level()));
  Index node_origin(node->X(), node->Y(), node->Z());
	
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
      auto clean_intersections = FilterClusters(start, direction, Normals, intersections);
      for (auto dist_isec: clean_intersections) {
	auto link_dist = dist_isec.first / norm;
	auto isec = dist_isec.second;
			  
	auto id = isec->second->id();
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


FluidTree::NodePtr SurfaceVoxeliser::ClassifyRegion(EdgeSiteTree::ConstNodePtr node) const {
  // Create output subtree
  auto ans = std::make_shared<FluidTree::Branch>(node->X(), node->Y(),
						 node->Z(), node->Level());
  const auto& dirs = Neighbours::GetDisplacements();

  // Iter over leaf nodes
  unsigned nLeafNodes = 0;
  node->IterDepthFirst(0, 0,
		       [&](EdgeSiteTree::ConstNodePtr edge_site_node) {
			 EdgeSitePtr edge_site = edge_site_node->Data();
			 auto& cuts = edge_site->closest_cut;
			 //assert(edge_site->Data())
			 // Examine all the links totting up how many are different types
			 unsigned nIn = 0;
			 unsigned nOut = 0;
			 unsigned nNA = 0;

			 FluidSite fsite;
			 for (auto i: range(26)) {
			   if (cuts[i].dist < 1) {
			     // there's a cut
			     auto& dir = dirs[i];
			     auto& faceId = cuts[i].id;
			     auto& norm = this->Normals[faceId];
			     auto dp = Vector::Dot(dir, norm);
			     if (dp > 0) {
			       // We are inside
			       ++nIn;

			       auto& label = this->Labels[faceId];
			       if (label < 0) {
				 // Wall
				 fsite.links[i].type = Intersection::Wall;
				 fsite.links[i].dist = cuts[i].dist;
			       } else {
				 // Iolet
				 // TODO: implement this properly
				 auto& iolet = this->Iolets[label];
				 
				 fsite.links[i].type = iolet.IsInlet ? Intersection::Inlet : Intersection::Outlet;
				 fsite.links[i].dist = cuts[i].dist;
				 fsite.links[i].id = iolet.Id;
			       }
			     } else {
			       ++nOut;
			       // This site better not have any inside!
			     }
			   } else {
			     // no cut
			     ++nNA;
			     fsite.links[i].type = Intersection::None;
			   }
			 }

			 // Now some sanity checks
			 assert((nIn + nOut + nNA) == 26);

			 if (nIn) {
			   // One or more interior sites, this one must be added to output
			   if (nOut) {
			     std::cout << *edge_site_node << std::endl;
			   }
			   assert (nOut == 0);

			   // Find closest wall link
			   float wall_dist = std::numeric_limits<float>::infinity();
			   int faceId;
			   for (auto i: range(26)) {
			     if (fsite.links[i].type == Intersection::Wall) {
			       if (fsite.links[i].dist < wall_dist) {
				 wall_dist = fsite.links[i].dist;
				 faceId = fsite.links[i].id;
			       }
			     }
			   }
			   if (wall_dist < 1) {
			     fsite.has_normal = true;
			     fsite.normal = this->Normals[faceId];
			   } else {
			     fsite.has_normal = false;
			   }

			   auto outsite = ans->GetCreate(edge_site_node->X(), edge_site_node->Y(), edge_site_node->Z(), 0);
			   outsite->Data().count = 1;
			   outsite->Data().leaf = std::make_shared<FluidSite>(fsite);
			   nLeafNodes += 1;
			 }
		       });
  ans->Data().count = nLeafNodes;
  return ans;
}
FluidTree SurfaceVoxeliser::operator()(const TriTree& inTree, const TriTree::Int tri_level) {
  typedef ParallelApply<FluidTree::NodePtr, TriTree::ConstNodePtr> Pool;

  Pool pool([this](TriTree::ConstNodePtr in) {
      auto edge_sites = ComputeIntersectionsForRegion(in);
      return ClassifyRegion(edge_sites);
    }, 1, std::numeric_limits<int>::max());

  std::deque<Pool::Future> result_queue;
  inTree.IterDepthFirst(tri_level, tri_level,
			[&](TriTree::ConstNodePtr node) {
			  result_queue.emplace_back(pool(node));
			});
  pool.Done();

  FluidTree outTree(inTree.Level());
	
  unsigned nLeafNodes = 0;
  while (!result_queue.empty()) {
    auto& fut = result_queue.front();

    auto edge_node = fut.get();
    nLeafNodes += edge_node->Data().count;
    outTree.Set(edge_node->X(), edge_node->Y(), edge_node->Z(), tri_level,
		edge_node);
    result_queue.pop_front();
  }
  outTree.Root()->Data().count = nLeafNodes;
  return outTree;
}
