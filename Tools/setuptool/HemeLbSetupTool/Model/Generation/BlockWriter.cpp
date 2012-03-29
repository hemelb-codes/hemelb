#include <cstdio>

#include <zlib.h>

#include "BlockWriter.h"
#include "GeometryWriter.h"
#include "Neighbours.h"
#include "GenerationError.h"

BlockWriter::BlockWriter(int blocksize) :
		maxBufferSize(geometry::GetMaxBlockRecordLength(blocksize)), buffer(
				new char[this->maxBufferSize]), writer(buffer, maxBufferSize) {
	this->Reset();
}

void BlockWriter::Reset() {
	nFluidSites = 0;
	CompressedBlockLength = 0;
	UncompressedBlockLength = 0;

	IsFinished = false;

	writer = hemelb::io::writers::xdr::XdrMemWriter(this->buffer,
			this->maxBufferSize);
}

BlockWriter::~BlockWriter() {
}

void BlockWriter::IncrementFluidSitesCount() {
	this->nFluidSites++;
}

void BlockWriter::Finish() {
	this->CompressedBlockLength = 0;
	this->UncompressedBlockLength = 0;

	if (this->nFluidSites > 0) {
		int ret; // zlib return code

		// How much data to compress?
		this->UncompressedBlockLength =
				this->writer.getCurrentStreamPosition();

		// Set up our compressor
		z_stream stream;
		stream.zalloc = Z_NULL;
		stream.zfree = Z_NULL;
		stream.opaque = Z_NULL;
		// Max compression
		ret = deflateInit(&stream, 9);
		if (ret != Z_OK)
			throw GenerationErrorMessage("Cannot init zlib structures");

		// Compute upper bound for how much space we'll need.
		int maxDeflatedLength = deflateBound(&stream,
				this->UncompressedBlockLength);
		char* compressedBuffer = new char[maxDeflatedLength];

		// Set input. The XDR buffer has to be char but zlib only works with
		// unsigned char. Just cast for now...
		stream.next_in = reinterpret_cast<unsigned char*>(this->buffer);
		stream.avail_in = this->UncompressedBlockLength;
		// Set output
		stream.next_out = reinterpret_cast<unsigned char*>(compressedBuffer);
		stream.avail_out = maxDeflatedLength;

		// Deflate. This should be it, if not their was an error.
		ret = deflate(&stream, Z_FINISH);
		if (ret != Z_STREAM_END)
			throw GenerationErrorMessage("Error compressing buffer");

		// How much space did we actually use?
		this->CompressedBlockLength = reinterpret_cast<char*>(stream.next_out)
				- compressedBuffer;

		// Tell zlib to clean up.
		ret = deflateEnd(&stream);
		if (ret != Z_OK)
			throw GenerationErrorMessage("Cannot free zlib structures");

		std::swap(buffer, compressedBuffer);
		delete[] compressedBuffer;
	} else {
		delete[] buffer;
		buffer = NULL;
	}

	IsFinished = true;
}

void BlockWriter::Write(GeometryWriter& gw) {
	if (this->nFluidSites > 0) {
		if (this->buffer == NULL)
			throw GenerationErrorMessage("Cannot write NULL buffer");

		std::fwrite(this->buffer, 1, this->CompressedBlockLength, gw.bodyFile);
	}
	*(gw.headerEncoder) << this->nFluidSites << this->CompressedBlockLength
			<< this->UncompressedBlockLength;
}
