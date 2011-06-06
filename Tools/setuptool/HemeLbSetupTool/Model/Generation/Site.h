#ifndef HEMELBSETUPTOOL_SITE_H
#define HEMELBSETUPTOOL_SITE_H

#include <vector>

#include "Index.h"

class Site {
public:
	unsigned int GetType() const;
	unsigned int GetConfig() const;
	Vector BoundaryNormal;
	double BoundaryDistance;
	bool IsEdge;
	Vector WallNormal;
	double WallDistance;
	std::vector<double> CutDistances;
protected:

};

#endif // HEMELBSETUPTOOL_SITE_H
