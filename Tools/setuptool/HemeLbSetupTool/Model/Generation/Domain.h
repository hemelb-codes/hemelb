#ifndef HEMELBSETUPTOOL_DOMAIN_H
#define HEMELBSETUPTOOL_DOMAIN_H

#include <vector>

#include "Index.h"
#include "GetSet.h"
#include "Site.h"
#include "BlockWriter.h"
class Block;
//class Site;
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

	GETTER(BlockSize, int);
	SETTER(BlockSize, int);

	GETTER(BlockCounts, Index);
	SETTER(BlockCounts, Index);

	GETTER(Origin, Vector);
	SETTER(Origin, Vector);

	GETTER(VoxelSize, double);
	SETTER(VoxelSize, double);

protected:
	Vector Origin;
	Index BlockCounts;
	unsigned int BlockSize;
	double VoxelSize;

	std::vector<Block*> blocks;

	inline Index* TranslateIndex(const int i) {
		Index* bInd = new Index();
		this->TranslateIndex(i, *bInd);
		return bInd;
	}
	inline void TranslateIndex(const int k, Index& ans) {
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

	inline int TranslateIndex(const Index& ijk) {
		return (ijk[0] * this->BlockCounts[1] + ijk[1]) * this->BlockCounts[2]
				+ ijk[2];
	}

	friend class BlockIterator;
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
