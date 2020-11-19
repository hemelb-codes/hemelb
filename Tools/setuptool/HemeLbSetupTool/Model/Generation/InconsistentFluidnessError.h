// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELBSETUPTOOL_INCONSISTENTFLUIDNESSERROR_H
#define HEMELBSETUPTOOL_INCONSISTENTFLUIDNESSERROR_H

#include <string>
#include "GenerationError.h"

class Site;

struct InconsistentFluidnessError: public GenerationError {
	InconsistentFluidnessError(const Site& s1, const Site& s2, const int nHits);
	~InconsistentFluidnessError() throw () {
	}
	virtual const char* what() const throw ();
	const Site& s1;
	const Site& s2;
	const int nHits;
private:
	std::string msg;
};

struct InconsistentIntersectRayError: public GenerationError {
  InconsistentIntersectRayError(const Site& s1, const Site& s2, const int nHits,
			     const bool sinside, const bool ninside);
	~InconsistentIntersectRayError() throw () {
	}
	virtual const char* what() const throw ();
	const Site& s1;
	const Site& s2;
	const int nHits;
	const bool sinside;
	const bool ninside;
private:
	std::string msg;
};

#endif // HEMELBSETUPTOOL_INCONSISTENTFLUIDNESSERROR_H
