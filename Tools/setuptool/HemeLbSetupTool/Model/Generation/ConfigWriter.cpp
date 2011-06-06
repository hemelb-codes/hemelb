#include "ConfigWriter.h"

#include "io/XdrFileWriter.h"
#include "io/XdrMemWriter.h"
#include "BlockWriter.h"

ConfigWriter::ConfigWriter(const std::string& OutputConfigFile, int StressType, int BlockSize,
			Index BlockCounts, double VoxelSize, Vector Origin) : OutputConfigFile(OutputConfigFile) {

	// Copy in key data
	//this->OutputConfigFile = OutputConfigFile;
	this->StressType = StressType;
	this->BlockSize = BlockSize;
	this->VoxelSize = VoxelSize;

	for (unsigned int i = 0; i < 3; ++i) {
		this->BlockCounts[i] = BlockCounts[i];
		this->Origin[i] = Origin[i];
	}

	// The encoder we're going to use for the body of the file. This needs to
	// be advanced to the right place in the file.
	this->bodyEncoder = new hemelb::io::XdrFileWriter(this->OutputConfigFile);
	// Do this by using this encoder to write the preceding parts; alias it to
	// a better name.
	hemelb::io::XdrWriter& headerEncoder = *this->bodyEncoder;

	// Preamble
	headerEncoder << this->StressType;

	// Blocks in each dimension
	for (unsigned int i = 0; i < 3; ++i)
		headerEncoder << this->BlockCounts[i];

	// Sites along 1 dimension of a block
	headerEncoder << this->BlockSize;

	// Voxel Size, in metres
	headerEncoder << this->VoxelSize;

	// Position of site index (0,0,0) in block index (0,0,0), in
	// metres in the STL file's coordinate system
	for (unsigned int i = 0; i < 3; ++i)
		headerEncoder << this->Origin[i];

	// (Dummy) Header

	// Note the start of the header
	this->headerStart = headerEncoder.getCurrentStreamPosition();
	unsigned int nBlocks = this->BlockCounts[0] * this->BlockCounts[1]
			* this->BlockCounts[2];

	// Write a dummy header
	for (unsigned int i = 0; i < nBlocks; ++i) {
		headerEncoder << 0 << 0;
	}

	this->bodyStart = headerEncoder.getCurrentStreamPosition();

	// Setup the encoder for the header
	this->headerBufferLength = this->bodyStart - this->headerStart;
	this->headerBuffer = new char[this->headerBufferLength];
	this->headerEncoder = new hemelb::io::XdrMemWriter(this->headerBuffer,
			this->headerBufferLength);

}

ConfigWriter::~ConfigWriter() {
	delete this->headerEncoder;
	delete this->headerBuffer;
	// Check this is still here as Close() will delete this
	if (this->bodyEncoder != NULL)
		delete this->bodyEncoder;
}

void ConfigWriter::Close() {
	// The hemelb::io::XdrFileWriter d'tor is currently the only way to close
	// the file
	delete this->bodyEncoder;
	this->bodyEncoder = NULL;

	// Reopen it, write the header buffer, close it.
	std::FILE* cfg = std::fopen(this->OutputConfigFile.c_str(), "r+");
	std::fseek(cfg, this->headerStart, SEEK_SET);
	std::fwrite(this->headerBuffer, 1, this->headerBufferLength, cfg);
	std::fclose(cfg);

}

BlockWriter* ConfigWriter::StartNextBlock() {
	return new BlockWriter(*this);
}

