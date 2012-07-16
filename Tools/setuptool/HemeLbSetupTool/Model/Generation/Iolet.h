// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

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
