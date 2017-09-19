//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELBSETUPTOOL_BUFFERPOOL_H
#define HEMELBSETUPTOOL_BUFFERPOOL_H
#include <stack>

// Allocates and frees or reuses buffers of a given size.
class BufferPool {
public:
	// C'tor- argument is the size of buffers to work with.
	BufferPool(unsigned int);
	~BufferPool();
	// Returns an uninitialized buffer
	char* New();
	// Return a buffer to the pool or free it
	void Free(char*);
	// Returns the size of buffers managed by the pool
	unsigned int GetSize() const;
private:
	unsigned int size;
	std::stack<char*> unused;
};

#endif
