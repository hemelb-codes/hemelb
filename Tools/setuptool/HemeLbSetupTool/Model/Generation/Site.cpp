#include <limits>

#include "Site.h"
#include "Block.h"
#include "Domain.h"

#include "Iolet.h"
#include "config.h"

#include "Neighbours.h"

// C'tor
Site::Site(Block& block, Index& index) :
	Position(block.GetDomain().CalcPositionFromIndex(index)), block(block),
			index(index) {
	this->Init();
}

// C'tor with index constructed in-place
Site::Site(Block& block, unsigned int i, unsigned int j, unsigned int k) :
	block(block), index(i, j, k) {
	this->Position = this->block.GetDomain().CalcPositionFromIndex(this->index);
	this->Init();
}

// private helper for common parts of c'tors
void Site::Init() {
	// Make sure this is false
	this->AdjacentIolet = NULL;
	// We need to find the closest.
	this->BoundaryDistance = std::numeric_limits<double>::infinity();
	this->WallDistance = std::numeric_limits<double>::infinity();

	this->CutDistances.resize(Neighbours::n,
			std::numeric_limits<double>::infinity());
	this->CutCellIds.resize(Neighbours::n, -1);
}

// Compute the type of the site from it's properties
unsigned int Site::GetType() const {
	if (this->IsFluid) {
		if (this->AdjacentIolet) {
			if (this->AdjacentIolet->IsInlet)
				return hemelb::INLET_TYPE;
			else
				return hemelb::OUTLET_TYPE;
		} else {
			return hemelb::FLUID_TYPE;
		}
	} else {
		return hemelb::SOLID_TYPE;
	}
}

// Compute the full config unsigned int
unsigned int Site::GetConfig() {
	unsigned int type = this->GetType();
	unsigned int cfg = type;
	// If solid, we're done
	if (!this->IsFluid)
		return cfg;

	// Fluid sites now
	if (type == hemelb::FLUID_TYPE && !this->IsEdge)
		// Simple fluid sites
		return cfg;

	// A complex one

	if (this->IsEdge) {
		// Bit fiddle the boundary config. See comment below
		// by BOUNDARY_CONFIG_MASK for definition.
		unsigned int boundary = 0;
		for (NeighbourIterator neighIt = this->beginall(); neighIt
				!= this->endall(); ++neighIt) {
			if (!neighIt->IsFluid) {
				// If the lattice vector is cut, set the flag
				boundary |= 1 << neighIt.GetNeighbourIndex();
			}
		}
		// Shift the boundary bit field to the appropriate
		// place and set these bits in the type
		cfg |= boundary << hemelb::BOUNDARY_CONFIG_SHIFT;

		if (boundary)
			// Set this bit if we've hit any solid sites
			cfg |= hemelb::PRESSURE_EDGE_MASK;

	}

	if (type != hemelb::FLUID_TYPE) {
		// It must be an inlet or outlet
		// Shift the index left and set the bits
		cfg |= this->BoundaryId << hemelb::BOUNDARY_ID_SHIFT;
	}

	return cfg;
}

// Get start and end LaterNeighbourIterators for this
LaterNeighbourIterator Site::begin() {
	return LaterNeighbourIterator(*this);
}
LaterNeighbourIterator Site::end() {
	return LaterNeighbourIterator(*this, Neighbours::nLater);
}
// Get start and end NeighbourIterators for this
NeighbourIterator Site::beginall() {
	return NeighbourIterator(*this);
}

NeighbourIterator Site::endall() {
	return NeighbourIterator(*this, Neighbours::n);
}

/*
 * Implement the NeighbourIterators
 */

NeighbourIteratorBase::NeighbourIteratorBase() {
}

// copy c'tor
NeighbourIteratorBase::NeighbourIteratorBase(const NeighbourIteratorBase& other) :
	site(other.site), domain(other.domain), i(other.i), maxI(other.maxI),
			index(other.index) {
}

void NeighbourIteratorBase::Init(Site& site, unsigned int startpos) {
	this->site = &site;
	this->index = site.index;
	this->domain = &site.block.domain;
	this->i = startpos;
	// Advance to the first extant neighbour
	while (this->i < this->maxI && !this->IsCurrentValid()) {
		++this->i;
	}
}

// C'tor for a given place, default is start
LaterNeighbourIterator::LaterNeighbourIterator(Site& site,
		unsigned int startpos) {
	this->maxI = Neighbours::nLater;
	this->Init(site, startpos);
}

// C'tor for a given place, default is start
NeighbourIterator::NeighbourIterator(Site& site, unsigned int startpos) {
	this->maxI = Neighbours::n;
	this->Init(site, startpos);
}

// assignment operator
NeighbourIteratorBase& NeighbourIteratorBase::operator=(
		const NeighbourIteratorBase& other) {
	if (this == &other) {
		return (*this);
	}
	this->site = other.site;
	this->domain = other.domain;
	this->i = other.i;
	this->maxI = other.maxI;

	return (*this);
}

// Is the current a site within the full domain?
bool NeighbourIteratorBase::IsCurrentValid() {
	this->index = this->site->index + this->GetVector();

	for (unsigned int j = 0; j < 3; ++j) {
		// If ind is out of the Domain, the current is invalid
		if (this->index[j] < 0 or this->index[j]
				>= this->domain->VoxelCounts[j])
			return false;
	}
	return true;
}

NeighbourIteratorBase& NeighbourIteratorBase::operator++() {
	// Note it is an error to increment an iterator past it's end, so we don't
	// need to handle that case.

	// Go to the next, then keep advancing if it's not valid
	++this->i;
	while (this->i < this->maxI && !this->IsCurrentValid()) {
		++this->i;
	}
	return *this;
}

// Test for equality with another
bool NeighbourIteratorBase::operator==(const NeighbourIteratorBase& other) const {
	return (other.site == this->site) && (other.i == this->i);
}
// Inequality
bool NeighbourIteratorBase::operator!=(const NeighbourIteratorBase& other) const {
	return !(*this == other);
}

// Dereference
NeighbourIteratorBase::reference NeighbourIteratorBase::operator*() {
	return this->domain->GetSite(this->index);
}
// Member lookup
NeighbourIteratorBase::pointer NeighbourIteratorBase::operator->() {
	return &(*(*this));
}

// Get the lattice vector for the current neighbour
Index LaterNeighbourIterator::GetVector() {
	return Neighbours::vectors[Neighbours::laterNeighbourIndices[this->i]];
}
// Get the lattice vector for the current neighbour
Index NeighbourIterator::GetVector() {
	return Neighbours::vectors[this->i];
}
