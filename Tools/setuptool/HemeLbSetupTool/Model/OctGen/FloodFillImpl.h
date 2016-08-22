// -*- mode: c++; -*-
#ifndef HEMELBSETUPTOOL_FLOODFILLIMPL_H
#define HEMELBSETUPTOOL_FLOODFILLIMPL_H
#include <array>
#include <exception>

#include "FluidSiteTree.h"

typedef std::array<FluidTree::Int, 4> Idx;
class StopIteration : public std::exception {
};

Idx GetStart(const FluidTree& tree);

class MaskBuilder;

// Alternative octree holding a binary mask and using arrays rather than
// pointers
class MaskTree {
public:
	typedef uint16_t Int;
	MaskTree(Int nl);
	bool Get(Int x, Int y, Int z);
private:
	static constexpr size_t NA() {
		return std::numeric_limits<size_t>::max();
	}

	friend MaskBuilder;
	friend class FloodFillTests;
	Int nLevels;
	std::vector<std::vector<size_t>> indices;
	std::vector<bool> mask;
};

// Class to build a MaskTree from a FluidTree
class MaskBuilder : public FluidTree::ConstVisitor {
public:
	typedef FluidTree::Int Int;

	MaskBuilder(const FluidTree& t);

	MaskTree operator()();

private:
	friend FluidTree::Node;
	friend class FloodFillTests;

	static size_t LocalIndex(FluidTree::ConstNodePtr np);
	virtual void Arrive(FluidTree::ConstNodePtr np);
	virtual void Depart(FluidTree::ConstNodePtr n) ;

    const FluidTree& tree;
    MaskTree ans;
    std::vector<size_t> cursors;
};

#endif
