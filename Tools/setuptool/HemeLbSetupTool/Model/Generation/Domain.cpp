#include "Site.h"
#include "Block.h"
#include "Domain.h"
#include "Debug.h"

Domain::Domain(double VoxelSizeMetres, double SurfaceBoundsWorking[6],
		unsigned int BlockSize) :
	BlockSize(BlockSize), VoxelSizeMetres(VoxelSizeMetres) {
	double min, max, size, extra, siteZero;
	int nSites, nBlocks, remainder, totalBlocks = 1;

	/*
	 * Here we are setting the location of our domain's origin in the input
	 * space and the number of sites along each axis. Sites will all have
	 * positions of:
	 * 		Origin + Index * VoxelSize,
	 * where:
	 * 		0 <= Index[i] < nSites[i]
	 *
	 * We also require that there be at least one solid site outside the fluid
	 * sites. For the case of axis-aligned faces which are an integer number
	 * of VoxelSizes apart (e.g. synthetic datasets!) this can cause numerical
	 * issues for the classifier if all the points that are "outside" are very
	 * close to the surface so we further require that these sites are a
	 * little further from the bounding box of the PolyData.
	 */
	for (unsigned int i = 0; i < 3; ++i) {
		// Bounds of the vtkPolyData
		min = SurfaceBoundsWorking[2 * i];
		max = SurfaceBoundsWorking[2 * i + 1];
		size = max - min;

		nSites = int(size / this->GetVoxelSizeWorking());
		/* Since int() truncates, we have:
		 * 		0 < size/VoxelSize - nSites < 1.
		 * Hence we need nSites + 1 links and therefore nSites + 2 sites
		 */
		nSites += 2;

		/* The extra distance from size to the distance from x[0] to x[nSites -1]
		 */
		extra = (nSites - 1) * this->GetVoxelSizeWorking() - size;

		/* To avoid numerical problems with the classifier, ensure that the
		 * sites just outside the fluid region are at least 1% of a VoxelSize
		 * away.
		 */
		if (extra < this->GetVoxelSizeWorking() / 100.) {
			// They weren't, so add one to the # sites and recalculate extra
			nSites += 1;
			extra = (nSites - 1) * this->GetVoxelSizeWorking() - size;
		}

		/* Now ensure this extra space is equally balanced before & after the
		 * fluid region with the placement of the first site.
		 */
		siteZero = min - 0.5 * extra;

		// Now work out how many blocks we require.
		nBlocks = nSites / BlockSize;
		remainder = nSites % BlockSize;
		if (remainder)
			++nBlocks;
		// Set the member vars for this axis
		this->OriginWorking[i] = siteZero;
		this->BlockCounts[i] = nBlocks;
		this->SiteCounts[i] = nBlocks * BlockSize;
		totalBlocks *= nBlocks;
	}
	// Resize the block vector
	this->blocks.resize(totalBlocks);
	Log() << "Domain size " << this->BlockCounts << std::endl;
}

Vector Domain::CalcPositionWorkingFromIndex(const Index& index) const {
	Vector ans(index);
	ans *= this->GetVoxelSizeWorking();
	ans += this->OriginWorking;
	return ans;
}

Block& Domain::GetBlock(const Index& index) {
	int i = TranslateIndex(index);
	Block* bp = blocks[i];
	// If the block hasn't been created yet, do so.
	if (bp == NULL) {
		bp = this->blocks[i] = new Block(*this, index, this->BlockSize);
	}
	return *bp;
}

Site& Domain::GetSite(const Index& gIndex) {
	Block& block = this->GetBlock(gIndex / this->BlockSize);
	return block.GetGlobalSite(gIndex);
}

BlockIterator Domain::begin() {
	return BlockIterator(*this);
}

BlockIterator Domain::end() {
	return BlockIterator(*this, Index(this->BlockCounts[0], 0, 0));
}

BlockIterator::BlockIterator() :
	domain(NULL), current(0, 0, 0) {
	this->maxima = this->domain->BlockCounts - 1;
}

BlockIterator::BlockIterator(Domain& dom) :
	domain(&dom), current(0, 0, 0) {
	this->maxima = this->domain->BlockCounts - 1;
}

BlockIterator::BlockIterator(Domain& dom, const Index& start) :
	domain(&dom), current(start) {
	this->maxima = this->domain->BlockCounts - 1;
}

BlockIterator::BlockIterator(const BlockIterator& other) :
	domain(other.domain), current(other.current), maxima(other.maxima) {
}
//	BlockIterator::~BlockIterator();

BlockIterator& BlockIterator::operator=(const BlockIterator& other) {
	if (this == &other) {
		return (*this);
	}
	this->domain = other.domain;
	this->current = other.current;
	this->maxima = other.maxima;

	return (*this);
}

BlockIterator& BlockIterator::operator++() {
	// Note it is an error to increment an iterator past it's end, so we don't
	// need to handle that case.
	int pos;
	// Delete any unnecessary blocks
	for (int i = this->current[0] - 1; i < this->current[0] + 1; ++i) {
		if (i < 0)
			continue;
		if (i == this->current[0] && i != this->maxima[0])
			continue;

		for (int j = this->current[1] - 1; j < this->current[1] + 1; ++j) {
			if (j < 0)
				continue;
			if (j == this->current[1] && j != this->maxima[1])
				continue;

			for (int k = this->current[2] - 1; k < this->current[2] + 1; ++k) {
				if (k < 0)
					continue;
				if (k == this->current[2] && k != this->maxima[2])
					continue;

				// This block can no longer be reached from the current or later
				// blocks, so delete, and set pointer to null
				pos = this->domain->TranslateIndex(i, j, k);
				delete this->domain->blocks[pos];
				this->domain->blocks[pos] = NULL;
			}
		}
	}

	// Update the index vector
	this->current[2] += 1;
	if (this->current[2] == this->domain->BlockCounts[2]) {
		this->current[2] = 0;

		this->current[1] += 1;
		if (this->current[1] == this->domain->BlockCounts[1]) {
			this->current[1] = 0;

			this->current[0] += 1;
		}
	}
	return *this;
}

bool BlockIterator::operator==(const BlockIterator& other) const {
	return (other.domain == this->domain) && (other.current == this->current);
}

bool BlockIterator::operator!=(const BlockIterator& other) const {
	return !(*this == other);
}

BlockIterator::reference BlockIterator::operator*() {
	return this->domain->GetBlock(this->current);
}

BlockIterator::pointer BlockIterator::operator->() {
	return &(*(*this));
}
