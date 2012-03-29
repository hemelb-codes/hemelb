#ifndef HEMELBSETUPTOOL_GEOMETRYWRITER_H
#define HEMELBSETUPTOOL_GEOMETRYWRITER_H

#include <string>
#include <cstdio>

#include "Index.h"
#include "BlockWriter.h"

#include "io/writers/xdr/XdrWriter.h"
using hemelb::io::writers::xdr::XdrWriter;

class GeometryWriter {
public:
	GeometryWriter(const std::string& OutputGeometryFile,
			int BlockSize, Index BlockCounts, double VoxelSizeMetres, Vector OriginMetres);

	~GeometryWriter();

	void Close();
	BlockWriter& StartNextBlock();

protected:
	std::string OutputGeometryFile;
	int BlockSize;
	Index BlockCounts;
	double VoxelSizeMetres;
	Vector OriginMetres;

	int headerStart;
	XdrWriter* headerEncoder;
	unsigned int headerBufferLength;
	char *headerBuffer;

	int bodyStart;
	FILE* bodyFile;
	BlockWriter blockWriter;

	friend class BlockWriter;
};

#endif // HEMELBSETUPTOOL_GEOMETRYWRITER_H
