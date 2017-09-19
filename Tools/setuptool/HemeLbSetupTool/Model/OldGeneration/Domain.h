// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELBSETUPTOOL_DOMAIN_H
#define HEMELBSETUPTOOL_DOMAIN_H

#include <vector>

#include "Index.h"
#include "GetSet.h"
#include "BlockWriter.h"
class Block;
class Site;
class BlockIterator;

// Note that the "working" units for this are voxels.

class Domain {
public:
	typedef BlockIterator iterator;
	/*
	 * C'tor
	 * SurfaceBounds - bounds of the surface, in standard VTK order
	 * (x_min, x_max, y_min, y_max, z_min, z_max), in voxels.
	 * BlockSize - number of sites along one dimension.
	 */
	Domain(double OriginWorking[3], unsigned SiteCounts[3], unsigned BlockSize = 8);

	Vector CalcPositionWorkingFromIndex(const Index& index) const;
	Block& GetBlock(const Index& index);
	Site& GetSite(const Index& index);
	BlockIterator begin();
	BlockIterator end();

	GETTER(BlockSize, int);SETTER(BlockSize, int);

	GETTER(BlockCounts, Index);SETTER(BlockCounts, Index);

	GETTER(SiteCounts, Index);

	/*
	 * These TranslateIndex member functions translate between 3d and 1a
	 * indices, i.e.
	 * 		ijk = (i * ny + j) * nz + k
	 * and the inverse. The direction is inferred from the type of the first
	 * argument.
	 */

	// 1d => 3d
	inline Index* TranslateIndex(const unsigned int i) {
		Index* bInd = new Index();
		this->TranslateIndex(i, *bInd);
		return bInd;
	}
	// 1d => 3d, putting the answer in an existing Index
	inline void TranslateIndex(const unsigned int k, Index& ans) {
		ans[2] = k % this->BlockCounts[2];
		int j = k / this->BlockCounts[2];

		ans[1] = j % this->BlockCounts[1];
		int i = j / this->BlockCounts[1];

		ans[0] = i % this->BlockCounts[0];
#ifdef CHECK_BOUNDS
		if (i / this->BlockCounts[0])
		throw IndexError;
#endif
	}

	// 3d => 1d
	inline int TranslateIndex(const Index& ijk) {
		return (ijk[0] * this->BlockCounts[1] + ijk[1]) * this->BlockCounts[2]
				+ ijk[2];
	}
	inline int TranslateIndex(const unsigned int& i, const unsigned int& j,
			const unsigned int& k) {
		return (i * this->BlockCounts[1] + j) * this->BlockCounts[2] + k;
	}

protected:
	Index BlockCounts;
	Index SiteCounts;
	Vector OriginWorking;
	unsigned int BlockSize;

	std::vector<Block*> blocks;

	friend class BlockIterator;
	friend class NeighbourIteratorBase;
};

class BlockIterator: public std::iterator<std::forward_iterator_tag, Block> {
public:
	BlockIterator();
	BlockIterator(Domain& dom);
	BlockIterator(Domain& dom, const Index& start);
	BlockIterator(const BlockIterator& other);
	//	~BlockIterator();

	BlockIterator& operator=(const BlockIterator& other);
	BlockIterator& operator++();
	bool operator==(const BlockIterator& other) const;
	bool operator!=(const BlockIterator& other) const;
	reference operator*();
	pointer operator->();

protected:
	Domain* domain;
	Index current;
	Index maxima;
};
#endif // HEMELBSETUPTOOL_DOMAIN_H
