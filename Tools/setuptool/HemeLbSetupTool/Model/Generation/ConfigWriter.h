#ifndef HEMELBSETUPTOOL_CONFIGWRITER_H
#define HEMELBSETUPTOOL_CONFIGWRITER_H

#include <string>
#include <cstdio>

#include "Index.h"

#include "io/writers/xdr/XdrWriter.h"
using hemelb::io::writers::xdr::XdrWriter;

class BlockWriter;

class ConfigWriter {
public:
	ConfigWriter(const std::string& OutputConfigFile, int StressType,
			int BlockSize, Index BlockCounts, double VoxelSize, Vector Origin);

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
	FILE* bodyFile;
	friend class BlockWriter;
};

#endif // HEMELBSETUPTOOL_CONFIGWRITER_H
