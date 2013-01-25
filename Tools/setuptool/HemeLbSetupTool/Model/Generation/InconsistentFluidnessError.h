//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

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

#endif // HEMELBSETUPTOOL_INCONSISTENTFLUIDNESSERROR_H
