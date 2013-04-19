// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELBSETUPTOOL_GEOMETRYWRITER_H
#define HEMELBSETUPTOOL_GEOMETRYWRITER_H

#include <string>
#include <cstdio>

#include "Index.h"

#include "io/writers/xdr/XdrWriter.h"
using hemelb::io::writers::xdr::XdrWriter;

class BlockWriter;
class BufferPool;

class GeometryWriter {
public:
	GeometryWriter(const std::string& OutputGeometryFile, int BlockSize,
			Index BlockCounts, double VoxelSizeMetres, Vector OriginMetres);

	~GeometryWriter();

	void Close();
	BlockWriter* StartNextBlock();

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
	BufferPool* BlockBufferPool;
	friend class BlockWriter;
};

#endif // HEMELBSETUPTOOL_GEOMETRYWRITER_H
