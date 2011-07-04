#include <cstdio>

#include "BlockWriter.h"
#include "ConfigWriter.h"

BlockWriter::BlockWriter(ConfigWriter &cfg) :
	configWriter(&cfg), nFluidSites(0) {
	this->maxBufferSize = (4 * 1 // unsigned int, site config
			+ 8 * 4 // doubles, boundary
			+ 8 * 4 // doubles, wall
			+ 8 * 14) // doubles, cut distances
			* cfg.BlockSize * cfg.BlockSize * cfg.BlockSize; // number of blocks
	this->buffer = new char[this->maxBufferSize];
	this->memWriter = new hemelb::io::XdrMemWriter(this->buffer,
			this->maxBufferSize);
}

BlockWriter::~BlockWriter() {
	delete this->memWriter;
	delete[] this->buffer;
}

void BlockWriter::IncrementFluidSitesCount() {
	this->nFluidSites++;
}

void BlockWriter::Finish() {
	unsigned int blockLength = 0;
	if (this->nFluidSites > 0) {
		// Work out how much space we've written to the buffer.
		blockLength = this->memWriter->getCurrentStreamPosition();

		// Write the buffer contents to the file.
		std::fwrite(this->buffer, 1, blockLength, this->configWriter->bodyFile);
	}

	// If there are no fluid sites, write nothing except two zeros to the header.
	// If there are fluid sites, blockLength was set above.
	(*this->configWriter->headerEncoder) << this->nFluidSites << blockLength;
}
