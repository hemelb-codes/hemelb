// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELBSETUPTOOL_BLOCK_H
#define HEMELBSETUPTOOL_BLOCK_H

#include <vector>

#include "Index.h"
class Site;
class vtkOBBTree;
#include "Domain.h"

typedef std::vector<Site*> SiteVec;
typedef SiteVec::iterator SiteIterator;
class Block {
public:
	//typedef std::vector<Site*>::iterator iterator;

	Block(Domain&, const Index&, const unsigned int&);
	~Block();

	Site& GetGlobalSite(const Index&);
	Site& GetLocalSite(const Index&);

	inline SiteIterator begin() {
		return this->sites.begin();
	}

	inline SiteIterator end() {
		return this->sites.end();
	}

	inline Domain& GetDomain() const {
		return this->domain;
	}
	inline const Index& GetIndex() const {
		return this->index;
	}
    vtkOBBTree * CreateOBBTreeModel(double extraSize) const;
    
    const Site & Middle() const {
        return *sites[sites.size()/2];
    }
protected:
	unsigned int size;
	const Index index;
	const Index min;
	const Index max;
	Domain& domain;
	SiteVec sites;

	inline unsigned int TranslateIndex(const Index& ind) {
		return (ind[0] * this->size + ind[1]) * this->size + ind[2];
	}
	friend class NeighbourIteratorBase;
	friend class LaterNeighbourIterator;
};

#endif // HEMELBSETUPTOOL_BLOCK_H
