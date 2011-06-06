#ifndef HEMELBSETUPTOOL_CONFIGWRITER_H
#define HEMELBSETUPTOOL_CONFIGWRITER_H

#include <string>

#include "Index.h"

namespace hemelb {
namespace io {
class XdrWriter;
}
}
using hemelb::io::XdrWriter;

class BlockWriter;

class ConfigWriter {
public:
	ConfigWriter(const std::string& OutputConfigFile, int StressType, int BlockSize,
			Index BlockCounts, double VoxelSize, Vector Origin);

	~ConfigWriter();

	void Close();
	BlockWriter* StartNextBlock();

protected:
	std::string OutputConfigFile;
	int StressType;
	int BlockSize;
	Index BlockCounts;
	double VoxelSize;
	Vector Origin;

	int headerStart;
	XdrWriter* headerEncoder;
	unsigned int headerBufferLength;
	char *headerBuffer;

	int bodyStart;
	XdrWriter* bodyEncoder;

	friend class BlockWriter;
};

#endif // HEMELBSETUPTOOL_CONFIGWRITER_H
