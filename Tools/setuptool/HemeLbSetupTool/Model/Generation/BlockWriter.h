#ifndef HEMELBSETUPTOOL_BLOCKWRITER_H
#define HEMELBSETUPTOOL_BLOCKWRITER_H

#include <stddef.h>
#include "ConfigWriter.h"
#include "io/writers/xdr/XdrMemWriter.h"

/*
 * Extension of a hemelb::io::XdrWriter that notes how many fluid sites, in how
 * much space, have been written. It then pushes this to the ConfigWriter's
 * headerEncoder and writes to the ConfigWriter's body
 */

class BlockWriter {
public:
	BlockWriter(ConfigWriter &cfg);
	~BlockWriter();

	void IncrementFluidSitesCount();

	void Finish();
	// Overload << to delegate to the XdrMemWriter
	template<typename T>
	BlockWriter& operator<<(T const & value) {
		(*this->memWriter) << value;
		return *this;
	}

	hemelb::io::writers::xdr::XdrMemWriter* memWriter;

protected:
	ConfigWriter* configWriter;
	unsigned int nFluidSites;
	size_t maxBufferSize;
	char *buffer;
};

#endif // HEMELBSETUPTOOL_BLOCKWRITER_H
