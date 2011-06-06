#ifndef HEMELBSETUPTOOL_GETSET_H
#define HEMELBSETUPTOOL_GETSET_H

// Macros for Get/Set methods
#define GETTER(name, type) inline type Get##name (void) {return this->name;}
#define SETTER(name, type) inline void Set##name(type val) {this->name = val;}

#endif // HEMELBSETUPTOOL_GETSET_H
