// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELBSETUPTOOL_GETSET_H
#define HEMELBSETUPTOOL_GETSET_H

// Macros for Get/Set methods
#define GETTER(name, type) inline type Get##name (void) {return this->name;}
#define SETTER(name, type) inline void Set##name(type val) {this->name = val;}

#endif // HEMELBSETUPTOOL_GETSET_H
