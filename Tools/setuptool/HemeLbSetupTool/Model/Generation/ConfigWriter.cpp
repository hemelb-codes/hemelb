#include "ConfigWriter.h"

#include "io/writers/xdr/XdrFileWriter.h"
#include "io/writers/xdr/XdrMemWriter.h"
#include "BlockWriter.h"

ConfigWriter::ConfigWriter(const std::string& OutputConfigFile, int StressType,
		int BlockSize, Index BlockCounts, double VoxelSize, Vector Origin) :
	OutputConfigFile(OutputConfigFile) {

	// Copy in key data
	//this->OutputConfigFile = OutputConfigFile;
	this->StressType = StressType;
	this->BlockSize = BlockSize;
	this->VoxelSize = VoxelSize;

	for (unsigned int i = 0; i < 3; ++i) {
		this->BlockCounts[i] = BlockCounts[i];
		this->Origin[i] = Origin[i];
	}

	{
		hemelb::io::writers::xdr::XdrFileWriter encoder(this->OutputConfigFile);

		// Preamble
		encoder << this->StressType;

		// Blocks in each dimension
		for (unsigned int i = 0; i < 3; ++i)
			encoder << this->BlockCounts[i];

		// Sites along 1 dimension of a block
		encoder << this->BlockSize;

		// Voxel Size, in metres
		encoder << this->VoxelSize;

		// Position of site index (0,0,0) in block index (0,0,0), in
		// metres in the STL file's coordinate system
		for (unsigned int i = 0; i < 3; ++i)
			encoder << this->Origin[i];

		// (Dummy) Header

		// Note the start of the header
		this->headerStart = encoder.getCurrentStreamPosition();
		unsigned int nBlocks = this->BlockCounts[0] * this->BlockCounts[1]
				* this->BlockCounts[2];

		// Write a dummy header
		for (unsigned int i = 0; i < nBlocks; ++i) {
			encoder << 0 << 0;
		}

		this->bodyStart = encoder.getCurrentStreamPosition();

		// Setup the encoder for the header
		this->headerBufferLength = this->bodyStart - this->headerStart;
		this->headerBuffer = new char[this->headerBufferLength];
		this->headerEncoder = new hemelb::io::writers::xdr::XdrMemWriter(this->headerBuffer,
				this->headerBufferLength);
	}
	// "encoder" declared in the block above will now have been destructed, thereby closing the file.
	// Reopen it for writing the blocks we're about to generate.
	this->bodyFile = std::fopen(this->OutputConfigFile.c_str(), "a");
}

ConfigWriter::~ConfigWriter() {
	delete this->headerEncoder;
	delete[] this->headerBuffer;
	// Check this is still here as Close() will delete this
	if (this->bodyFile != NULL)
		std::fclose(this->bodyFile);
}

void ConfigWriter::Close() {
	// Close the config
	std::fclose(this->bodyFile);
	this->bodyFile= NULL;

	// Reopen it, write the header buffer, close it.
	std::FILE* cfg = std::fopen(this->OutputConfigFile.c_str(), "r+");
	std::fseek(cfg, this->headerStart, SEEK_SET);
	std::fwrite(this->headerBuffer, 1, this->headerBufferLength, cfg);
	std::fclose(cfg);

}

BlockWriter* ConfigWriter::StartNextBlock() {
	return new BlockWriter(*this);
}

