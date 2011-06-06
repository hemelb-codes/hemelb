#include "BlockWriter.h"
#include "ConfigWriter.h"

BlockWriter::BlockWriter(ConfigWriter &cfg) :
	configWriter(&cfg) {
	this->count = 0;
	this->blockStart
			= this->configWriter->bodyEncoder->getCurrentStreamPosition();
}

void BlockWriter::IncrementFluidSitesCount() {
	this->count++;
}

void BlockWriter::Finish() {
	// Work out how much space our block occupies
	unsigned int blockLength =
			this->configWriter->bodyEncoder->getCurrentStreamPosition()
					- this->blockStart;

	// Write this to the header
	(*this->configWriter->headerEncoder) << this->count << blockLength;

}
