// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELBSETUPTOOL_GETSET_H
#define HEMELBSETUPTOOL_GETSET_H

// Macros for Get/Set methods
#define GETTER(name, type) inline type Get##name (void) {return this->name;}
#define SETTER(name, type) inline void Set##name(type val) {this->name = val;}

#endif // HEMELBSETUPTOOL_GETSET_H
