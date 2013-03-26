// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELBSETUPTOOL_CYLINDERGENERATOR_H
#define HEMELBSETUPTOOL_CYLINDERGENERATOR_H

#include "GeometryGenerator.h"

#include "Index.h"
#include "GetSet.h"
#include "Iolet.h"
#include "GenerationError.h"

class GeometryWriter;
class Site;
class BlockWriter;

struct CylinderData {
	// Cylinder parameters
	Vector Centre;
	Vector Axis;
	double Radius;
	double Length;
};

class CylinderGenerator : public GeometryGenerator {
public:
	CylinderGenerator();
	virtual ~CylinderGenerator();

	inline void SetCylinderCentre(Vector v) {
		this->Cylinder->Centre = v;
	}
	inline void SetCylinderAxis(Vector n) {
		this->Cylinder->Axis = n;
	}
	inline void SetCylinderRadius(double r) {
		this->Cylinder->Radius = r;
	}
	inline void SetCylinderLength(double l) {
		this->Cylinder->Length = l;
	}

private:
	virtual void ComputeBounds(double []) const;
	void ClassifySite(Site& site);
	void ComputeCylinderNormalAtAPoint(Vector& wallNormal, const Vector& surfacePoint, const Vector& cylinderAxis) const;
	CylinderData* Cylinder;
protected:
	virtual int BlockInsideOrOutsideSurface(const Block &block) {
        return 0;
	}
};

#endif // HEMELBSETUPTOOL_CYLINDERGENERATOR_H
