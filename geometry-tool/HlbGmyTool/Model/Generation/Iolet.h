// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELBSETUPTOOL_IOLET_H
#define HEMELBSETUPTOOL_IOLET_H

#include "Index.h"
struct Iolet {
  Vector Centre;
  Vector Normal;

  double Radius;
  int Id;
  bool IsInlet;
};

#endif  // HEMELBSETUPTOOL_IOLET_H
