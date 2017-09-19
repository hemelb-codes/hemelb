// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELBSETUPTOOL_SQUAREDUCTGENERATOR_H
#define HEMELBSETUPTOOL_SQUAREDUCTGENERATOR_H

#include "GeometryGenerator.h"

#include "Index.h"
#include "GetSet.h"
#include "Iolet.h"
#include "GenerationError.h"

class GeometryWriter;
class Site;
class BlockWriter;

struct SquareDuctData {
	Vector LowerBound;
	Vector UpperBound;
	int OpenAxis;
};

class SquareDuctGenerator : public GeometryGenerator {
public:
	SquareDuctGenerator();
	virtual ~SquareDuctGenerator();

	inline void SetLowerBound(Vector v) {
		this->SquareDuct->LowerBound = v;
	}
	inline void SetUpperBound(Vector v) {
			this->SquareDuct->UpperBound = v;
		}
	inline void SetOpenAxis(int i) {
		this->SquareDuct->OpenAxis = i;
	}

private:
	virtual void ComputeBounds(double []) const;
	void ClassifySite(Site& site);
	SquareDuctData* SquareDuct;
protected:
	virtual int BlockInsideOrOutsideSurface(const Block &block) {
        return 0;
    }
};

#endif // HEMELBSETUPTOOL_SQUAREDUCTGENERATOR_H
