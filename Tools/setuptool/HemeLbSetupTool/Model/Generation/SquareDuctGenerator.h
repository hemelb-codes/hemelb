// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

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
