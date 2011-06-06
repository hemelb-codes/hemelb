#ifndef HEMELBSETUPTOOL_BLOCKWRITER_H
#define HEMELBSETUPTOOL_BLOCKWRITER_H

#include "io/XdrWriter.h"
using hemelb::io::XdrWriter;

#include "ConfigWriter.h"

/*
 * Extension of a hemelb::io::XdrWriter that delegates all it's writing to
 * the ConfigWriter's bodyEncoder, but notes how many fluid sites, in how
 *  much space, have been written. It then pushes this to the ConfigWriter's
 *  headerEncoder.
 */

class BlockWriter : public XdrWriter {
public:
	BlockWriter(ConfigWriter &cfg);
	void IncrementFluidSitesCount();

	void Finish();

protected:
	ConfigWriter* configWriter;
	unsigned int count;
	unsigned int blockStart;

    // Methods to write basic types to the Xdr object.
	template <typename T>
	void _write(T& val) {
		this->configWriter->bodyEncoder->_write(val);
	}
};


#endif // HEMELBSETUPTOOL_BLOCKWRITER_H
