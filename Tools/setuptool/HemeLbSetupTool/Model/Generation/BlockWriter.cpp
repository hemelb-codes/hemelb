#include <cstdio>

#include <zlib.h>

#include "BlockWriter.h"
#include "GeometryWriter.h"
#include "Neighbours.h"

BlockWriter::BlockWriter(GeometryWriter &cfg) :
	geometryWriter(&cfg), nFluidSites(0) {
	this->maxBufferSize = geometry::GetMaximumBlockRecordLength(cfg.BlockSize);
	this->buffer = new char[this->maxBufferSize];
	this->memWriter = new hemelb::io::writers::xdr::XdrMemWriter(this->buffer,
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
	unsigned int compressedBlockLength = 0;
	unsigned int uncompressedBlockLength = 0;
	if (this->nFluidSites > 0) {
		int ret; // zlib return code

		// How much data to compress?
		uncompressedBlockLength = this->memWriter->getCurrentStreamPosition();

		// Set up our compressor
		z_stream stream;
		stream.zalloc = Z_NULL;
		stream.zfree = Z_NULL;
		stream.opaque = Z_NULL;
	    // Max compression
		ret = deflateInit(&stream, 9);
		// assert ret == Z_OK

		// Compute upper bound for how much space we'll need.
		int maxDeflatedLength = deflateBound(&stream, uncompressedBlockLength);
		unsigned char* compressedBuffer = new unsigned char[maxDeflatedLength];

		// Set input. The XDR buffer has to be char but zlib only works with
		// unsigned char. Just cast for now...
		stream.next_in = reinterpret_cast<unsigned char*>(this->buffer);
		stream.avail_in = uncompressedBlockLength;
		// Set output
		stream.next_out = compressedBuffer;
		stream.avail_out = maxDeflatedLength;

		// Deflate. This should be it, if not their was an error.
		ret = deflate(&stream, Z_FINISH);
		// assert ret == Z_STREAM_END

		// How much space did we actually use?
		compressedBlockLength = stream.next_out - compressedBuffer;

		// Tell zlib to clean up.
		ret = deflateEnd(&stream);
		// assert ret == Z_OK

		// Write the buffer contents to the file.
		std::fwrite(compressedBuffer, 1, compressedBlockLength,
				this->geometryWriter->bodyFile);
		delete[] compressedBuffer;
	}

	// If there are no fluid sites, write nothing except two zeros to the header.
	// If there are fluid sites, blockLength was set above.
	(*this->geometryWriter->headerEncoder) << this->nFluidSites
			<< compressedBlockLength << uncompressedBlockLength;
}
