#ifndef HEMELBSETUPTOOL_DOMAIN_H
#define HEMELBSETUPTOOL_DOMAIN_H

#include <vector>

#include "Index.h"
#include "GetSet.h"
//#include "Site.h"
#include "BlockWriter.h"
class Block;
class Site;
class BlockIterator;

class Domain {
public:
	typedef BlockIterator iterator;
	/*
	 * C'tor
	 * VoxelSize - voxel size, in metres
	 * SurfaceBounds - bounds of the surface, in standard VTK order
	 * (x_min, x_max, y_min, y_max, z_min, z_max), in metres.
	 * BlockSize - number of sites along one dimension.
	 */
	Domain(double VoxelSize, double SurfaceBounds[6],
			unsigned int BlockSize = 8);

	Vector CalcPositionFromIndex(const Index& index) const;
	Block& GetBlock(const Index& index);
	Site& GetSite(const Index& index);
	BlockIterator begin();
	BlockIterator end();

	GETTER(BlockSize, int);SETTER(BlockSize, int);

	GETTER(BlockCounts, Index);SETTER(BlockCounts, Index);

	GETTER(Origin, Vector);SETTER(Origin, Vector);

	GETTER(VoxelSize, double);SETTER(VoxelSize, double);

protected:
	Vector Origin;
	Index BlockCounts;
	Index VoxelCounts;
	unsigned int BlockSize;
	double VoxelSize;

	std::vector<Block*> blocks;
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
		ans.z = k % this->BlockCounts.z;
		int j = k / this->BlockCounts.z;

		ans.y = j % this->BlockCounts.y;
		int i = j / this->BlockCounts.y;

		ans.x = i % this->BlockCounts.x;
#ifdef CHECK_BOUNDS
		if (i / this->BlockCounts[0])
		throw IndexError;
#endif
	}

	// 3d => 1d
	inline int TranslateIndex(const Index& ijk) {
		return (ijk.x * this->BlockCounts.y + ijk.y) * this->BlockCounts.z
				+ ijk.z;
	}
	inline int TranslateIndex(const unsigned int& i, const unsigned int& j,
			const unsigned int& k) {
		return (i * this->BlockCounts.y + j) * this->BlockCounts.z + k;
	}

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
