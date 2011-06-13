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
	const Index GetDomainBlockCount();
	const int GetDomainBlockSize();

protected:
	void Init();
	Block& block;
	Index index;
	friend class NeighbourIteratorBase;
	friend class LaterNeighbourIterator;
};

// Base for iterating over neighbours
class NeighbourIteratorBase: public std::iterator<std::forward_iterator_tag,
		Site> {
public:
	NeighbourIteratorBase(Site& site, unsigned int startpos = 0);
	NeighbourIteratorBase(const NeighbourIteratorBase& other);

	NeighbourIteratorBase& operator=(const NeighbourIteratorBase& other);
	NeighbourIteratorBase& operator++();
	bool operator==(const NeighbourIteratorBase& other) const;
	bool operator!=(const NeighbourIteratorBase& other) const;
	reference operator*();
	pointer operator->();

	virtual unsigned int GetNeighbourIndex();

protected:
	Site* site;
	Domain* domain;
	unsigned int i;

	void AdvanceToValid();
	bool IsCurrentInDomain();
	virtual bool IsCurrentValid() = 0;
	Index GetVector();
};

// Iterator for getting all the later (i.e. further on in memory) neighbouring sites of a given site
class LaterNeighbourIterator: public NeighbourIteratorBase {
public:
	inline LaterNeighbourIterator(Site& site, unsigned int startpos = 0) :
		NeighbourIteratorBase(site, startpos) {
		this->AdvanceToValid();
	}
protected:
	bool IsCurrentValid();
};

// Iterator for getting all the later (i.e. further on in memory) neighbouring sites of a given site
class NeighbourIterator: public NeighbourIteratorBase {
public:
	inline NeighbourIterator(Site& site, unsigned int startpos = 0) :
		NeighbourIteratorBase(site, startpos) {
		this->AdvanceToValid();
	}
protected:
	bool IsCurrentValid();
};
#endif // HEMELBSETUPTOOL_SITE_H
