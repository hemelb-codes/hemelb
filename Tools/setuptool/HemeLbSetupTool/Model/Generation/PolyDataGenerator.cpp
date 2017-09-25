#include "PolyDataGenerator.h"
#include "vtkPolyData.h"

std::string PolyDataGenerator::GetOutputGeometryFile(void) {
  return this->OutputGeometryFile;
}
void PolyDataGenerator::SetOutputGeometryFile(std::string val) {
  this->OutputGeometryFile = val;
}

std::vector<Iolet*>& PolyDataGenerator::GetIolets() {
  return this->Iolets;
}
void PolyDataGenerator::SetIolets(std::vector<Iolet*> iv) {
  this->Iolets = std::vector<Iolet*>(iv);
}

void PolyDataGenerator::SetOriginWorking(double x, double y, double z) {
  this->OriginWorking[0] = x;
  this->OriginWorking[1] = y;
  this->OriginWorking[2] = z;
}

void PolyDataGenerator::SetSiteCounts(unsigned x, unsigned y, unsigned z) {
  this->SiteCounts[0] = x;
  this->SiteCounts[1] = y;
  this->SiteCounts[2] = z;
}

void PolyDataGenerator::GetSeedPointWorking(double out[3]) {
  for (unsigned int i = 0; i < 3; ++i)
    out[i] = this->SeedPointWorking[i];
  return;
}
void PolyDataGenerator::SetSeedPointWorking(double out[3]) {
  for (unsigned int i = 0; i < 3; ++i)
    this->SeedPointWorking[i] = out[i];
}
void PolyDataGenerator::SetSeedPointWorking(double x, double y, double z) {
  this->SeedPointWorking[0] = x;
  this->SeedPointWorking[1] = y;
  this->SeedPointWorking[2] = z;
}

vtkPolyData* PolyDataGenerator::GetClippedSurface(void) {
  return this->ClippedSurface;
}
void PolyDataGenerator::SetClippedSurface(vtkPolyData* val) {
  this->ClippedSurface = val;
}

void PolyDataGenerator::Execute() throw (GenerationError) {
  std::cout << "Implement Execute!!" << std::endl;
}
