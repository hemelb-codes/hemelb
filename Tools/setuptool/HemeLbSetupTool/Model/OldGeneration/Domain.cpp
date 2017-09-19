//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#include "Site.h"
#include "Block.h"
#include "Domain.h"
#include "Debug.h"

Domain::Domain(double OriginWorking[3], unsigned SiteCounts[3], unsigned BlockSize) :
	BlockSize(BlockSize) {

	int remainder, totalBlocks = 1;

	for (unsigned int i = 0; i < 3; ++i) {
		// Copy in
		this->OriginWorking[i] = OriginWorking[i];
		this->SiteCounts[i] = SiteCounts[i];

		// Now work out how many blocks we require.
		this->BlockCounts[i] = this->SiteCounts[i] / BlockSize;
		remainder = this->SiteCounts[i] % BlockSize;
		if (remainder)
			++this->BlockCounts[i];

		// Adjust the site count
		this->SiteCounts[i] = this->BlockCounts[i] * BlockSize;

		totalBlocks *= this->BlockCounts[i];
	}
	// Resize the block vector
	this->blocks.resize(totalBlocks);
	Log() << "Domain size " << this->BlockCounts << std::endl;
}

Vector Domain::CalcPositionWorkingFromIndex(const Index& index) const {
	Vector ans(index);
	ans += this->OriginWorking;
	return ans;
}

Block& Domain::GetBlock(const Index& index) {
	int i = this->TranslateIndex(index);
	Block* bp = this->blocks[i];
	// If the block hasn't been created yet, do so.
	if (!bp) {
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
