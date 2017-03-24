#include <exception>
#include <deque>

#include "FloodFill.h"

#include "util/Vector3D.h"
#include "Neighbours.h"
#include "range.hpp"

#define unpack(idx) idx[0], idx[1], idx[2], idx[3]

Idx GetStart(const FluidTree& tree) {
  Idx ans;
  try {
    tree.IterDepthFirst(0,0,
			[&ans](FluidTree::ConstNodePtr n) {
			  ans = {n->X(), n->Y(), n->Z(), 0};
			  throw StopIteration();
			});
  } catch (StopIteration& stop) {
    return ans;
  }
  throw std::runtime_error("Could not get a start point from tree");
}


#define PrintIdx(i) std::cout << "[" << i[0] << ", " << i[1] << ", "<< i[2] << ", "<< i[3] << "]" << std::endl

MaskTree FloodFill(const FluidTree& tree) {
  const auto& dirs = Neighbours::GetDisplacements();
  typedef std::deque<Idx> Queue;
  
  auto seed = GetStart(tree);
    
  MaskTree seen(tree.Level());
  seen.GetCreate(unpack(seed));
  
  // Stack holds sites that are def fluid but we don't know about their
  // neighbours
  Queue stack;
  stack.push_back(seed);
  auto BULK = std::make_shared<FluidSite>();
  
  while (!stack.empty()) {
    // Grab the point from the queue
    auto pt = stack.back();
    stack.pop_back();
    
    auto leaf_node_ptr = tree.Get(unpack(pt));
    if (leaf_node_ptr) {
      // Node exists already: it must be fluid and have been created
      // by SurfaceVoxeliser. Look at all links and queue all sites
      // with intersection == none
      const auto& node_data = leaf_node_ptr->Data();
      const auto& links = node_data.leaf->links;
      for (auto i: range(26)) {
	if (links[i].type == Intersection::None) {
	  // enqueue (pt + dirs[i]) if we have not seen it
	  Idx neigh{FluidTree::Int(pt[0]+dirs[i][0]),
	      FluidTree::Int(pt[1]+dirs[i][1]),
	      FluidTree::Int(pt[2]+dirs[i][2]),
	      pt[3]};

	  if (!seen.Get(unpack(neigh))) {
	    PrintIdx(neigh);
	    seen.GetCreate(unpack(neigh));
	    stack.emplace_back(neigh);
	  }
	}
      }
    } else {
      // Node does not exist, but it must be fluid (as we are only
      // enqueuing nodes that are linked to a fluid site without cuts)

      // Could create a leaf node for it with the basic fluid site
      // state, but that is a lot of work

      // Instead, just iter over all the neighbours and enqueue if not
      // seen
      for (auto i: range(26)) {
	Idx neigh{FluidTree::Int(pt[0]+dirs[i][0]),
	    FluidTree::Int(pt[1]+dirs[i][1]),
	    FluidTree::Int(pt[2]+dirs[i][2]),
	    pt[3]};
	if (!seen.Get(unpack(neigh))) {
	  seen.GetCreate(unpack(neigh));
	  stack.emplace_back(neigh);
	}
      }
      
    }
  }
  return seen;
}




