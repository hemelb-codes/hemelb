// -*- mode: c++; -*-
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

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

  void SetIolets(std::vector<Iolet*> iv);

  void SetOriginWorking(double x, double y, double z);

  void SetNumberOfLevels(int n);
  void SetTriangleLevel(int n);
  
  void GetSeedPointWorking(double out[3]);
  void SetSeedPointWorking(double out[3]);
  void SetSeedPointWorking(double x, double y, double z);

  vtkPolyData* GetClippedSurface(void);
  void SetClippedSurface(vtkPolyData* val);

private:
  // Members set from outside to initialise
  double OriginWorking[3];
  uint16_t NumberOfLevels;
  uint16_t TriangleLevel;
  std::string OutputGeometryFile;
  std::vector<Iolet> Iolets;
  double SeedPointWorking[3];
  vtkPolyData* ClippedSurface;
  
};

#endif
