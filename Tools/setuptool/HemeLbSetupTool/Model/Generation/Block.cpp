#include "Index.h"

#include "Site.h"
#include "Block.h"
#include "Domain.h"

Block::Block(Domain& dom, const Index& ind, const unsigned int& size) :
	size(size), index(ind), min(ind * size), max((ind + 1) * size), domain(dom) {
	this->sites.resize(size * size * size);
	unsigned int ijk = 0;
	for (unsigned int i = ind[0] * size; i < (ind[0] + 1) * size; ++i) {
		for (unsigned int j = ind[1] * size; j < (ind[1] + 1) * size; ++j) {
			for (unsigned int k = ind[2] * size; k < (ind[2] + 1) * size; ++k) {
				this->sites[ijk] = new Site(*this, i, j, k);
				++ijk;
			}
		}
	}
}

Block::~Block() {
	SiteIterator end = this->sites.end();
	SiteIterator current = this->sites.begin();
	for (; current != end; ++current) {
		// current will dereference to a Site*
		delete *current;
	}
}

Site& Block::GetGlobalSite(const Index& globalInd) {
	bool local = true;
	for (unsigned int i = 0; i < 3; ++i) {
		if (globalInd[i] < this->min[i] || globalInd[i] >= this->max[i]) {
			local = false;
			break;
		}
	}
	if (local)
		return this->GetLocalSite(globalInd - this->min);

	// Check if the coords belong to another block, i.e. any of
	// the local ones outside the range [0, self.size)

	return this->domain.GetSite(globalInd);
}

Site& Block::GetLocalSite(const Index& ind) {
	/*
	 * Get the site, creating it if it didn't exist.
	 */
	unsigned int ijk = this->TranslateIndex(ind);
	//	if (this->sites[ijk] == NULL)
	//		this->sites[ijk] = new Site(*this);
	return *this->sites[ijk];
}
