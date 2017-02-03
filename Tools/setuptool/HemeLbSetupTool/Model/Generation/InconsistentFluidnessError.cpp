//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#include "InconsistentFluidnessError.h"
#include "Site.h"
#include <sstream>

void FormatSite(std::ostringstream& msg, const Site& site) {
	msg << "site index " << site.GetIndex() << ", position " << site.Position
			<< ", which is ";
	if (site.IsFluid)
		msg << "fluid";
	else
		msg << "solid";
}


void FormatBoolInside(std::ostringstream& msg, const bool sinside) {
	if (sinside == true)
		msg << "fluid";
	else
		msg << "solid";
}

InconsistentFluidnessError::InconsistentFluidnessError(const Site& s1,
		const Site& s2, const int nHits) :
	s1(s1), s2(s2), nHits(nHits) {
	std::ostringstream msgStream;
	msgStream << "Inconsistent fluidness detected between ";
	FormatSite(msgStream, s1);
	msgStream << " and ";
	FormatSite(msgStream, s2);
	msgStream << " but found " << nHits << " intersections with the surface.";
	this->msg = msgStream.str();
}

InconsistentIntersectRayError::InconsistentIntersectRayError(const Site& s1,
	   const Site& s2, const int nHits, const bool sinside, const bool ninside) :
	   s1(s1), s2(s2), nHits(nHits), sinside(sinside), ninside(ninside) {
	std::ostringstream msgStream;
	msgStream << "Inconsistent fluidness from ray and intersect between ";
	FormatSite(msgStream, s1);
	msgStream << " and ";
	FormatSite(msgStream, s2);
	msgStream << " found " << nHits << 
	  " intersections with the surface, but according to Raycast site 1 is "; 
	FormatBoolInside(msgStream, sinside);
	msgStream  << " while site 2 is ";
	FormatBoolInside(msgStream, ninside);
	this->msg = msgStream.str();
}

const char* InconsistentFluidnessError::what() const throw () {
	return this->msg.c_str();
}


const char* InconsistentIntersectRayError::what() const throw () {
	return this->msg.c_str();
}
