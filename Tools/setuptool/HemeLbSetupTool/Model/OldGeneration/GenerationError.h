// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELBSETUPTOOL_GENERATIONERROR_H
#define HEMELBSETUPTOOL_GENERATIONERROR_H
#include <exception>
#include <string>

struct GenerationError: public std::exception {
	virtual const char* what() const throw () {
		return "GenerationError";
	}
};

struct GenerationErrorMessage: public GenerationError {
	GenerationErrorMessage(const std::string errorMessage) :
		msg(errorMessage) {
	}
	~GenerationErrorMessage() throw () {
	}
	virtual const char* what() const throw () {
		return msg.c_str();
	}

	const std::string msg;
};

#endif
