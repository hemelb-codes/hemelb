#include "Block.h"
#include "Index.h"
#include "Domain.h"
#include "Site.h"

Block::Block(Domain& dom, const Index& ind, unsigned int& size) :
	size(size), index(ind), min(ind * size), max((ind + 1) * size), domain(dom) {
	this->sites.resize(size * size * size);
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
	return this->sites[this->TranslateIndex(ind)];
}
