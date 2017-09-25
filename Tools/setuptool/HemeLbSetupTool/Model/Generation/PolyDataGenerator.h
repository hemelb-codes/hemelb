// -*- mode: c++; -*-
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELBSETUPTOOL_POLYDATAGENERATOR_H
#define HEMELBSETUPTOOL_POLYDATAGENERATOR_H

#include <vector>
#include "GenerationError.h"

class vtkPolyData;
class Iolet;

class PolyDataGenerator {
public:
  PolyDataGenerator();
  
  void Execute() throw (GenerationError);
  
  std::string GetOutputGeometryFile(void);
  void SetOutputGeometryFile(std::string val);

  std::vector<Iolet*>& GetIolets();
  void SetIolets(std::vector<Iolet*> iv);

  void SetOriginWorking(double x, double y, double z);

  void SetSiteCounts(unsigned x, unsigned y, unsigned z);

  void GetSeedPointWorking(double out[3]);
  void SetSeedPointWorking(double out[3]);
  void SetSeedPointWorking(double x, double y, double z);

  vtkPolyData* GetClippedSurface(void);
  void SetClippedSurface(vtkPolyData* val);

private:
  // Members set from outside to initialise
  double OriginWorking[3];
  unsigned SiteCounts[3];
  std::string OutputGeometryFile;
  std::vector<Iolet*> Iolets;
  double SeedPointWorking[3];
  vtkPolyData* ClippedSurface;
  
};

#endif
