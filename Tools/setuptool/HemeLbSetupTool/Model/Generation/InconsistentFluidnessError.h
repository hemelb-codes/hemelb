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
