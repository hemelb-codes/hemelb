#ifndef HEMELBSETUPTOOL_BLOCKWRITER_H
#define HEMELBSETUPTOOL_BLOCKWRITER_H

#include "ConfigWriter.h"

/*
 * Extension of a hemelb::io::XdrWriter that delegates all it's writing to
 * the ConfigWriter's bodyEncoder, but notes how many fluid sites, in how
 * much space, have been written. It then pushes this to the ConfigWriter's
 * headerEncoder.
 */

class BlockWriter {
public:
	BlockWriter(ConfigWriter &cfg);
	void IncrementFluidSitesCount();

	void Finish();
	// Overload << to delegate to the config writer.
	template<typename T>
	BlockWriter& operator<<(T const & value) {
		(*this->configWriter->bodyEncoder) << value;
		return *this;
	}

protected:
	ConfigWriter* configWriter;
	unsigned int count;
	unsigned int blockStart;

};

#endif // HEMELBSETUPTOOL_BLOCKWRITER_H
