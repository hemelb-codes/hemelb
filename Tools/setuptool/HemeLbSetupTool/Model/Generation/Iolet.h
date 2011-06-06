#ifndef HEMELBSETUPTOOL_IOLET_H
#define HEMELBSETUPTOOL_IOLET_H

//#include <vector>
#include "Index.h"
struct Iolet {
	//  std::vector<double> Centre;
	//  std::vector<double> Normal;
	Vector Centre;
	Vector Normal;

	double Radius;
	int Id;
	bool IsInlet;
};

//typedef std::vector<Iolet> IoletVector;

#endif // HEMELBSETUPTOOL_IOLET_H
