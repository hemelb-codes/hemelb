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


const char* InconsistentFluidnessError::what() const throw () {
	return this->msg.c_str();
}
