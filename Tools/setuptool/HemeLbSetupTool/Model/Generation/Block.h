#ifndef HEMELBSETUPTOOL_BLOCK_H
#define HEMELBSETUPTOOL_BLOCK_H

#include <vector>

#include "Index.h"
//#include "Site.h"
class Site;
class Domain;

class Block {
public:
	typedef std::vector<Site>::iterator iterator;

	Block(Domain&, const Index&, unsigned int&);
	Site& GetGlobalSite(const Index&);
	Site& GetLocalSite(const Index&);

	inline iterator begin() {
		return this->sites.begin();
	}

	inline iterator end() {
		return this->sites.end();
	}

protected:
	unsigned int size;
	const Index& index;
	const Index min;
	const Index max;
	Domain& domain;
	std::vector<Site> sites;

	inline unsigned int TranslateIndex(const Index& ind) {
		return (ind[0] * this->size + ind[1]) * this->size + ind[2];
	}
};

#endif // HEMELBSETUPTOOL_BLOCK_H
