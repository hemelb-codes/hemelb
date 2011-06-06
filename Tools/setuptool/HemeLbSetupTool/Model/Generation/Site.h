#ifndef HEMELBSETUPTOOL_SITE_H
#define HEMELBSETUPTOOL_SITE_H

#include <vector>

#include "Index.h"
class Iolet;
class Block;
class Domain;

class LaterNeighbourIterator;
class NeighbourIterator;

// A single lattice site
class Site {
public:
	Site(Block& block, Index& index);
	Site(Block& block, unsigned int i, unsigned int j, unsigned int k);
	unsigned int GetType() const;
	unsigned int GetConfig();
	bool IsFluidKnown;
	bool IsFluid;
	bool IsEdge;
	Vector BoundaryNormal;
	double BoundaryDistance;
	Vector WallNormal;
	double WallDistance;
	std::vector<double> CutDistances;
	std::vector<int> CutCellIds;
	Iolet* AdjacentIolet;
	unsigned int BoundaryId;
	Vector Position;

	LaterNeighbourIterator begin();
	LaterNeighbourIterator end();
	NeighbourIterator beginall();
	NeighbourIterator endall();
	inline const Index& GetIndex() {
		return this->index;
	}
protected:
	void Init();
	Block& block;
	Index index;
	friend class NeighbourIteratorBase;
};

// Base for iterating over neighbours
class NeighbourIteratorBase: public std::iterator<std::forward_iterator_tag,
		Site> {
public:
	NeighbourIteratorBase();
	NeighbourIteratorBase(const NeighbourIteratorBase& other);

	NeighbourIteratorBase& operator=(const NeighbourIteratorBase& other);
	NeighbourIteratorBase& operator++();
	bool operator==(const NeighbourIteratorBase& other) const;
	bool operator!=(const NeighbourIteratorBase& other) const;
	reference operator*();
	pointer operator->();

	virtual unsigned int GetNeighbourIndex() = 0;

protected:
	Site* site;
	Domain* domain;
	unsigned int i;
	unsigned int maxI;
	Index index;

	// Does the work of the constructor
	void Init(Site& site, unsigned int startpos);
	bool IsCurrentValid();
	virtual Index GetVector() = 0;
};

// Iterator for getting all the later (i.e. further on in memory) neighbouring sites of a given site
class LaterNeighbourIterator: public NeighbourIteratorBase {
public:
	LaterNeighbourIterator(Site& site, unsigned int startpos = 0);
	unsigned int GetNeighbourIndex();
protected:
	Index GetVector();
};

// Iterator for getting all the later (i.e. further on in memory) neighbouring sites of a given site
class NeighbourIterator: public NeighbourIteratorBase {
public:
	NeighbourIterator(Site& site, unsigned int startpos = 0);
	unsigned int GetNeighbourIndex();
protected:
	Index GetVector();
};
#endif // HEMELBSETUPTOOL_SITE_H
