#ifndef HEMELBSETUPTOOL_CONFIGGENERATOR_H
#define HEMELBSETUPTOOL_CONFIGGENERATOR_H

#include <string>
#include <vector>

// VTK bits we need
class vtkPolyData;
class vtkOBBTree;
class vtkPoints;
class vtkIdList;

#include "GetSet.h"
#include "Iolet.h"

class ConfigWriter;
class Site;

class ConfigGenerator {
public:
	ConfigGenerator();
	~ConfigGenerator();
	void Execute();
	void ClassifySite(Site& site);
	bool GetIsFluid(Site& site);

	inline double GetVoxelSize(void) {
		return this->VoxelSize;
	}
	inline void SetVoxelSize(double val) {
		this->VoxelSize = val;
	}

	inline std::string GetOutputConfigFile(void) {
		return this->OutputConfigFile;
	}
	inline void SetOutputConfigFile(std::string val) {
		this->OutputConfigFile = val;
	}

	inline std::vector<Iolet*>& GetIolets() {
		return this->Iolets;
	}
	inline void SetIolets(std::vector<Iolet*> iv) {
		this->Iolets = std::vector<Iolet*>(iv);
	}

	inline void GetSeedPoint(double out[3]) {
		for (unsigned int i = 0; i < 3; ++i)
			out[i] = this->SeedPoint[i];
		return;
	}
	inline void SetSeedPoint(double out[3]) {
		for (unsigned int i = 0; i < 3; ++i)
			this->SeedPoint[i] = out[i];
	}
	inline void SetSeedPoint(double x, double y, double z) {
		this->SeedPoint[0] = x;
		this->SeedPoint[1] = y;
		this->SeedPoint[2] = z;
	}

	inline vtkPolyData* GetClippedSurface(void) {
		return this->ClippedSurface;
	}
	inline void SetClippedSurface(vtkPolyData* val) {
		this->ClippedSurface = val;
	}

protected:
	// Members set from outside to initialise
	double VoxelSize;
	std::string OutputConfigFile;
	std::vector<Iolet*> Iolets;
	double SeedPoint[3];
	vtkPolyData* ClippedSurface;
	vtkOBBTree* Locator;

	// Members used internally
	vtkPoints* hitPoints;
	vtkIdList* hitCellIds;
	bool IsFirstSite;
};

#endif // HEMELBSETUPTOOL_CONFIGGENERATOR_H
