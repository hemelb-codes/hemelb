#include <exception>
#include <deque>

#include "FloodFill.h"
#include "FloodFillImpl.h"

#include "util/Vector3D.h"
#include "Neighbours.h"

//typedef hemelb::util::Vector3D<FluidTree::Int> Idx;

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

size_t LIndex(MaskTree::Int x, MaskTree::Int y, MaskTree::Int z,
		MaskTree::Int level) {
	MaskTree::Int lx = (x >> level) & 1;
	MaskTree::Int ly = (y >> level) & 1;
	MaskTree::Int lz = (z >> level) & 1;
	return (lx << 2) + (ly << 1) + lz;
}


MaskTree::MaskTree(Int nl) : nLevels(nl),  indices(nl), mask() {
	indices[nl-1] = std::vector<size_t>(8, NA());
}

bool MaskTree::Get(Int x, Int y, Int z) {
	auto cur_level = nLevels;
	size_t offset = 0;
	while (--cur_level) {
		Int lInd = LIndex(x, y, z, cur_level);
		offset = indices[cur_level][offset + lInd];
		if (offset == NA()) {
			return false;
		}
	}
	Int lInd = LIndex(x, y, z, cur_level);
	return mask[offset + lInd];
}

MaskBuilder::MaskBuilder(const FluidTree& t) : tree(t), ans(t.Level()), cursors(t.Level(), 0) {
}

MaskTree MaskBuilder::operator()() {
	tree.Root()->Accept(*this);
	return ans;
}

size_t MaskBuilder::LocalIndex(FluidTree::ConstNodePtr np) {
	return LIndex(np->X(), np->Y(), np->Z(), np->Level());
}

void MaskBuilder::Arrive(FluidTree::ConstNodePtr np) {
	Int level = np->Level();
	if (level == tree.Level()) {
		// Top level node implicily exists and points to offset zero of the
		// level below. So do nothing.
	} else if (level > 0) {
		// Non-leaf node
		auto i = cursors[level];
		size_t n;
		if (level == 1)
			n = ans.mask.size();
		else
			n = ans.indices[level - 1].size();

		ans.indices[level][i + LocalIndex(np)] = n;
		cursors[level - 1] = n;
		if (level == 1)
			ans.mask.resize(n + 8, false);
		else
			ans.indices[level - 1].resize(n + 8, MaskTree::NA());
	} else {
		// Must be leaf node
		auto i = cursors[level];
		ans.mask[i + LocalIndex(np)] = true;
	}
}

void MaskBuilder::Depart(FluidTree::ConstNodePtr n) {}


void FloodFill(FluidTree& tree) {
	Idx seed = GetStart(tree);
	const auto& dirs = Neighbours::GetDisplacements();

	typedef std::deque<Idx> Queue;

	// Stack holds sites that are def fluid but we don't know about their
	// neighbours
	Queue stack;
	stack.push_back(seed);

	while (!stack.empty()) {
		auto pt = stack.back();
		stack.pop_back();

		auto site = tree.Get(unpack(pt));

	}

}


