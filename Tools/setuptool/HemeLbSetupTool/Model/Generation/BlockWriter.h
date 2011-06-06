#ifndef HEMELBSETUPTOOL_BLOCKWRITER_H
#define HEMELBSETUPTOOL_BLOCKWRITER_H

#include "io/XdrWriter.h"
using hemelb::io::XdrWriter;

#include "ConfigWriter.h"

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
