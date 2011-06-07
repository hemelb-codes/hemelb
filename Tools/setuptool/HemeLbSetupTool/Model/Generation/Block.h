#ifndef HEMELBSETUPTOOL_BLOCK_H
#define HEMELBSETUPTOOL_BLOCK_H

#include <vector>

#include "Index.h"
//#include "Site.h"
class Site;
class Domain;

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

	inline Domain& GetDomain() {
		return this->domain;
	}
	inline const Index& GetIndex() {
		return this->index;
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
};

#endif // HEMELBSETUPTOOL_BLOCK_H
